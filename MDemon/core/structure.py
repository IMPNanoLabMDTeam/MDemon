import functools
import numbers
from copy import copy

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from .. import _PARTICLES, _STRUCTURE_NAMES, _STRUCTURES, _TOPOLOGIES
from ..lib.util import flat


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


class _StructureAttrContainer:
    _SETATTR_WHITELIST = ["atom_graph"]

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
            attrname = attr.import_attrname(sid)

            # getter and setter is defined in a loop,
            # so we need make_getter and make_setter
            # to pass them correct sid.
            def make_getter(sid):
                def getter(self):
                    return attr.__getitem__(self, sid)

                return getter

            def make_setter(sid):
                def setter(self, values):
                    return attr.__setitem__(self, sid, values)

                return setter

            setattr(
                cls, attrname, property(make_getter(sid), make_setter(sid), None, None)
            )

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


class _ImmutableBase:
    """
    Class used to shortcut :meth:`__new__` to :meth:`object.__new__`.

    When mixed via _StructureAttrContainer._mix this class has MRO priority. \
    Setting __new__ like this will avoid having to go through the \
    cache lookup if the class is reused.
    """

    __new__ = object.__new__


class _MutableBase:
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
    def family(self):
        return self._u.families[self._fname]

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

    def atom_connection(self):
        sid = ("Atom_" + self._fname, "Atom_" + self._fname)
        source = self._u._database.connection._source_register[sid]
        return source.values

    def atoms2graph(self):
        if "atom_graph" not in self.__dict__:
            atoms = self._u.atoms

            mtrx_connect = self.atom_connection()

            if isinstance(self, Atom):
                atm_list = self.ix
            else:
                atm_list = flat(self.atms)
            atm_set = set(atm_list)
            pairs = []
            for i in atm_list:
                for j in atoms[i].neighbors:
                    if j in atm_set and i > j:
                        pairs.append((int(i), int(j), {"label": mtrx_connect[i, j]}))

            G = nx.Graph()
            G.add_edges_from(pairs)
            self.atom_graph = G

    def draw_atoms(self, style, k=1):
        atoms = self._u.atoms

        if style == "graph":
            self.atoms2graph()
            G = self.atom_graph

            class_colors = {1: "red", 2: "green", 3: "blue", 4: "purple"}

            num_nodes = len(G.nodes)
            plt.figure(figsize=(num_nodes // 10, num_nodes // 10), dpi=300)
            pos = nx.spring_layout(G, k=k * np.sqrt(1 / num_nodes))

            node_colors = [class_colors[atoms[node].species] for node in G.nodes]

            nx.draw_networkx_nodes(
                G, pos, node_color=node_colors, node_size=250, alpha=0.8
            )
            nx.draw_networkx_edges(G, pos, alpha=0.5)
            nx.draw_networkx_labels(G, pos, font_size=10)

            edge_labels = nx.get_edge_attributes(G, "label")
            nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=5)

            plt.title(f"{self.sname} {self.ix}")

    def create_rings(self, multiring=False):
        """
        Create rings from an atom graph.
        """
        self.atoms2graph()

        cycles = list(nx.simple_cycles(self.atom_graph, length_bound=10))
        H = nx.Graph()
        # 添加简单环中的节点和边到子图 H 中
        for cycle in cycles:
            # 每个环中的节点
            H.add_nodes_from(cycle)
            # 每个环中的边
            H.add_edges_from(
                [(cycle[i], cycle[i + 1]) for i in range(len(cycle) - 1)]
                + [(cycle[-1], cycle[0])]
            )

        cycles = list(nx.minimum_cycle_basis(H))

        sname = "Ring_" + self._fname
        sname_atm = "Atom_" + self._fname
        sname_bnd = "Bond_" + self._fname

        # 获取每个回路的边的 label 属性
        def get_cycle_edge_labels(G, cycles):
            cycle_edge_labels = []
            for cycle in cycles:
                labels = []
                for i in range(len(cycle)):
                    u = cycle[i]
                    v = cycle[(i + 1) % len(cycle)]

                    labels.append(G[u][v]["label"])
                cycle_edge_labels.append(labels)
            return cycle_edge_labels

        # 获取回路的边的 label 属性, which is ix of bonds.
        cycle_edge_labels = get_cycle_edge_labels(self.atom_graph, cycles)

        self._u._database.register_structure(
            sname, {sname_atm: cycles, sname_bnd: cycle_edge_labels}
        )

        self._u._database.composition.to_connection(sname)

        if multiring:
            self._u._database.connection.union(
                "Ring_Base", "Multiring_Base", process_attr=True
            )

    @property
    def n_atoms(self):
        return len(self.atms)


class Particle(Structure, metaclass=_ParticleMeta):
    abbreviation = "pr"

    def prune(self):
        atoms = self.family.atoms

        n_effective_neighbors_dic = {}
        for ix in self.atms:
            i = 0
            a = atoms[ix]
            for ix1 in a.neighbors:
                if ix1 in self.atms:
                    i += 1
            n_effective_neighbors_dic[ix] = i

        visited_atms = []
        wait_to_visit_atms = copy(self.atms)
        finish = False
        while not finish:
            finish = True
            for ix in wait_to_visit_atms:
                if ix not in visited_atms:
                    if n_effective_neighbors_dic[ix] <= 1:
                        a = atoms[ix]
                        visited_atms.append(ix)
                        for ix1 in a.neighbors:
                            n_effective_neighbors_dic[ix1] -= 1
                        finish = False
        ixs = []
        for ix in self.atms:
            if ix not in visited_atms:
                ixs.append(ix)

        return Atom(np.array(ixs, dtype=np.int32), self.universe)


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


class Multiring(Topology):
    abbreviation = "mrng"

    @property
    def n_rings(self):
        return len(self.rngs)
