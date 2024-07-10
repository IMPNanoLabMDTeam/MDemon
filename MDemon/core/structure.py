import functools
import numbers

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from .. import _PARTICLES, _STRUCTURE_NAMES, _STRUCTURES, _TOPOLOGIES
from ..lib.util import asiterable


def flat(nums):
    """
    Flat the multi-layer nesting list.
    """
    res = []
    nums = asiterable(nums)
    for i in nums:
        if isinstance(i, list):
            res.extend(flat(i))
        else:
            res.append(i)
    return res


def _check_family_consistency(func):
    """
    In most of the time, `structure` objects \
    only interact with other objects/classes from \
    their own family, as a result we should check \
    whether those objects belong to the same family.
    """

    @functools.wraps(func)
    def wrapped(strucs):
        families = []
        for s in strucs:
            families.append(s.family)
        if len(set(families)) != 1:
            raise ValueError("Can't operate on objects from different families！")
        return func(strucs)

    return wrapped


class _ParticleMeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        _STRUCTURE_NAMES[name] = cls
        _PARTICLES.append(cls)
        _STRUCTURES.append(cls)


class _TopologyMeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        _STRUCTURE_NAMES[name] = cls
        _TOPOLOGIES.append(cls)
        _STRUCTURES.append(cls)


class _StructureAttrContainer(object):
    _SETATTR_WHITELIST = ["pairs"]

    @classmethod
    def _subclass(cls):
        newcls = type(cls.__name__, (cls,), {})

        return newcls

    @classmethod
    def _mix(cls, other, family):
        name = other.__name__
        fname = family.name
        newcls = type(name, (_ImmutableBase, cls, other), {})
        newcls._fname = fname
        newcls.family = family
        return newcls

    @classmethod
    def _add_prop(cls, attr):
        sids = attr.import_sids(cls)

        for sid in sids:

            def getter(self):
                return attr.__getitem__(self, sid)

            def setter(self, values):
                return attr.__setitem__(self, sid, values)

            attrname = attr.import_attrname(sid)
            setattr(cls, attrname, property(getter, setter, None, None))

    def __setattr__(self, attr, value):
        # `ag.this = 42` calls setattr(ag, 'this', 42)
        if not (
            attr.startswith("_")  # 'private' allowed
            or attr in self._SETATTR_WHITELIST  # known attributes allowed
            or hasattr(self, attr)  # preexisting (eg properties) allowed
        ):
            raise AttributeError(f"Cannot set arbitrary attributes to a {type(self)}")
        # if it is, we allow the setattr to proceed by deferring to the super
        # behaviour (ie do it)
        super(_StructureAttrContainer, self).__setattr__(attr, value)


class _ImmutableBase(object):
    """
    Class used to shortcut :meth:`__new__` to :meth:`object.__new__`.

    When mixed via _StructureAttrContainer._mix this class has MRO priority. \
    Setting __new__ like this will avoid having to go through the \
    cache lookup if the class is reused.
    """

    __new__ = object.__new__


class _MutableBase(object):
    """
    Base class that merges appropriate :class:`_StructureAttrContainer` classes.\n
    In it the instantiating class is fetched from :attr:`Universe._classes`.\
    The classes themselves are used as the cache dictionary keys for\
    simplicity in cache retrieval.
    """

    def __new__(cls, *args, **kwargs):
        # This pre-initialization wrapper must be pretty generic to
        # allow for different initialization schemes of the possible classes.
        # All we really need here is to fish a universe out of the arg list.
        try:
            u = args[-1].universe
        except (IndexError, AttributeError):
            try:
                # older AtomGroup init method..
                u = args[0][0].universe
            except (TypeError, IndexError, AttributeError):
                errmsg = (
                    f"No universe, or universe-containing object "
                    f"passed to the initialization of {cls.__name__}"
                )
                raise TypeError(errmsg) from None
        try:
            _fname = kwargs["family"]
        except KeyError:
            _fname = "Base"
        try:
            return object.__new__(u.families[_fname]._classes[cls])
        except KeyError:
            errmsg = f"class {cls.__name__} is not defined in the family {_fname}"
            raise TypeError(errmsg) from None


class Structure(_MutableBase):
    ## Rule of Abbreviation ########
    ################################
    # 1. Use lowercase letters.
    # 2. (The first consonant of decoration part) + The first main syllable without vowels.
    # For example: mrng for Multiring and strc for Structure.
    # 3. If the abbreviation is too short, plus the second main syllable without vowels.
    # For example: tp for Topology.
    # 4. If the fisrt main syllable has no consonant, combine the fisrt vowel
    # with the second main syllable without vowels.
    # For example: agl for angle.
    # 5. The definition of main syllable: the combination of "consonants + vowels +
    # (consonants without following vowels)" which has basic meanings of the structure.
    # For example: ring in multiring.
    abbreviation = "strc"
    _compo = []  # "composition"
    _fname = "Base"  # family name
    detailed = True  # control the verbosity of output

    def __init__(self, *args, **kwargs):
        try:
            if len(args) == 1:
                # List of structures.
                # Make sure they belong to same family!
                self._compo = args[0]
                self._compo2ix()
                u = self.compo[0].universe
            elif len(args) == 2:
                # ix : :attr:`index`
                ix, u = args
                self._ix = np.asarray(flat(ix), dtype=np.intp)
        except (
            AttributeError,  # couldn't find ix/universe
            TypeError,
        ):  # couldn't iterate the object we got
            errmsg = (
                "Can only initialise a Structure from an iterable of Particle/"
                "Topology objects eg: Atom([Molecule1, Cluster2, Ring3]) "
                "or an iterable of indices and a Universe reference "
                "eg: Atom([6], u)."
            )
            raise TypeError(errmsg) from None

        # indices for the objects I hold
        self._u = u

    def __getitem__(self, item):
        # supports
        # - integer access
        # - boolean slicing
        # - fancy indexing
        # because our _ix attribute is a numpy array
        # it can be sliced by all of these already,
        # so just return ourselves sliced by the item
        if item is None:
            raise TypeError("None cannot be used to index a group.")
        elif isinstance(item, numbers.Integral):
            return self.__class__(self.ix[item], self.universe)
        else:
            if isinstance(item, list) and item:  # check for empty list
                # hack to make lists into numpy arrays
                # important for boolean slicing
                item = np.array(item)
            # We specify _derived_class instead of self.__class__ to allow
            # subclasses, such as UpdatingAtomGroup, to control the class
            # resulting from slicing.
            return self.__class__(self.ix[item], self.universe)

    @_check_family_consistency
    def ix2ix(self, newcls):
        """
        Convert :attr:`index` to indices of \
        another :class:`Structure`.

        Parameters
        ----------
        newcls : `Structure`
            Targeted subclass of :class:`Structure`

        Returns
        -------
        indices : `numpy.ndarray`
        """
        abbrev = newcls.abbreviation
        try:
            return self.__getattribute__(abbrev + "s")
        except AttributeError:
            return self.__getattribute__(abbrev + "_ix")

    def _compo2ix(self):
        """
        Derive the indices from :attr:`_compo`.
        """
        ix = []
        for s in self._compo:
            ix.extend(flat(s.ix2ix(self.abbreviation)))
        self._ix = np.asarray(list(set(ix)), dtype=np.intp)

    @property
    def universe(self):
        """The underlying :class:`~MaxwellDemon.core.universe.Universe` \
        the structure belongs to.
        """
        return self._u

    @property
    def ix(self):
        """Unique indices of object's components."""
        if self._u.silent:
            return self._ix[self.silence]
        else:
            return self._ix

    @property
    def sname(self):
        """
        Structural name of this class.
        """
        return self.__class__.__name__ + "_" + self._fname

    @classmethod
    def import_sname(cls):
        """
        Why `class property` will not be supported \
        after Python 3.13?
        """
        return cls.__name__ + "_" + cls._fname

    def __len__(self):
        return len(self._ix)

    def atoms2graph(self):
        atoms = self._u.atoms
        if "pairs" not in self.__dict__:
            self.pairs = []
            for i in self.atms:
                for j in atoms[i].neighbors:
                    if j in self.atms and i > j:
                        self.pairs.append((i, j))
        return self.pairs

    def draw_atoms(self, style):
        atoms = self._u.atoms

        if style == "graph":
            pairs = self.atoms2graph()

            G = nx.Graph()
            G.add_edges_from(pairs)

            class_colors = {1: "red", 2: "green", 3: "blue", 4: "purple"}

            num_nodes = len(G.nodes)
            plt.figure(figsize=(min(20, num_nodes // 2), min(20, num_nodes // 2)))
            pos = nx.spring_layout(G, k=np.sqrt(1 / num_nodes))

            node_colors = [class_colors[atoms[node].species] for node in G.nodes]

            nx.draw_networkx_nodes(
                G, pos, node_color=node_colors, node_size=500, alpha=0.8
            )
            nx.draw_networkx_edges(G, pos, alpha=0.5)
            nx.draw_networkx_labels(G, pos)

            plt.title(f"{self.sname} {self.ix}")


class Particle(Structure, metaclass=_ParticleMeta):
    abbreviation = "pr"


class Atom(Particle):
    abbreviation = "atm"


class Group(Particle):
    abbreviation = "grp"


class Molecule(Particle):
    abbreviation = "mle"


class Topology(Structure, metaclass=_TopologyMeta):
    abbreviation = "tp"


class Chain(Topology):
    abbreviation = "chn"


class Bond(Chain):
    abbreviation = "bnd"


class Angle(Chain):
    abbreviation = "agl"


class Dihedral(Chain):
    abbreviation = "dh"


class Tree(Topology):
    abbreviation = "tr"


class Improper(Tree):
    abbreviation = "ipr"


class Ring(Topology):
    abbreviation = "rng"
