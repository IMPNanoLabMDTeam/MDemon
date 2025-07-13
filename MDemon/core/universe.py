import numpy as np

from .. import _PARTICLES, _TOPOLOGIES, _UNIVERSE_ATTRS
from .rw import get_reader_for
from .structure import Structure, _StructureAttrContainer


def _make_bases():
    bases = {}
    # The 'SBase' middle man is needed so that a single structureattr
    # patching applies automatically to all structures.
    SBase = bases[Structure] = _StructureAttrContainer._subclass()
    PBase = SBase._subclass()
    for cls in _PARTICLES:
        bases[cls] = PBase._subclass()
    TBase = SBase._subclass()
    for cls in _TOPOLOGIES:
        bases[cls] = TBase._subclass()

    return bases


def _make_classes(family):
    # Initializes the class cache.
    u = family._u
    classes = {}
    for cls in u._class_bases:
        classes[cls] = u._class_bases[cls]._mix(cls, family)
    family._classes = classes


def _generate_from_database(dbase):
    """
    Generate family-specific Universe versions of each :class:`Structure`,
    including :class:`Atom`, :class:`Bond`, :class:`Molecule` etc.
    """
    family = dbase.base
    _make_classes(family)

    # Put Structure level stuff from database into class
    # and attach attrs to universe.
    for attr in dbase.attrs:
        if attr.__class__ not in _UNIVERSE_ATTRS:
            family._process_attr(attr)
        else:
            dbase._u._add_prop(attr)

    # Generate literally everthing.
    family.instancing()


def _database_from_file_like(*inputfiles, **kwargs):
    for file_ in inputfiles:
        reader = get_reader_for(file_)
        with reader(file_) as r:
            database = r.parse(**kwargs)
            kwargs["database"] = database

    return database


class Universe:
    def __init__(self, *inputfiles, **kwargs) -> None:
        self.timestep = 0  # initial timestep
        self._s = False  # silent
        self._class_bases = _make_bases()
        self._database = _database_from_file_like(*inputfiles, **kwargs)
        self._database._u = self

        _generate_from_database(self._database)

    @property
    def silent(self):
        return self._s

    @silent.setter
    def silent(self, is_silent):
        if not isinstance(is_silent, bool):
            raise TypeError("The value of silence should be True or False.")
        self._s = is_silent

    @classmethod
    def _add_prop(cls, attr):
        """
        In fact, universe is just a bigger structure.\
        In current version, we don't consider multiple\
        universes, so only `self.timestep` is needed.
        """

        def getter(self):
            return attr.__getitem__(self.timestep)

        def setter(self, values):
            return attr.__setitem__(self.timestep, values)

        setattr(cls, attr.name, property(getter, setter, None, None))

    @property
    def universe(self):
        return self

    @property
    def families(self):
        return self._database.families


class UniverseAttrMeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        _UNIVERSE_ATTRS.append(cls)


class UniverseAttr:
    def __init__(self, *values, database):
        self.values = np.array(values)
        self._database = database
        database.add_Attr(self)

    def __getitem__(self, tix):
        try:
            return self.values[tix]
        except IndexError:
            return self.values[0]


class Timestep(UniverseAttr, metaclass=UniverseAttrMeta):
    name = "timestep"


class Box(UniverseAttr, metaclass=UniverseAttrMeta):
    name = "box"
