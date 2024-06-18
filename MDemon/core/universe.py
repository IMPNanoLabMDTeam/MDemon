from .. import _PARTICLES, _TOPOLOGIES
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


def _generate_from_database(family):
    """Generate family-specific Universe versions of each :class:`Structure`,
    including :class:`Atom`, :class:`Bond`, :class:`Molecule` etc.
    """
    _make_classes(family)

    # Put Structure level stuff from database into class
    for attr in family._u._database.attrs:
        family._process_attr(attr)

    # Generate literally everthing.
    family.instancing()


def _database_from_file_like(*inputfiles, **kwargs):
    file0 = inputfiles[0]
    reader = get_reader_for(file0)
    with reader(file0) as r:
        database = r.parse(**kwargs)

    return database


class Universe(object):
    def __init__(self, *inputfiles) -> None:
        self._s = False  # silent
        self._class_bases = _make_bases()
        self._database = _database_from_file_like(*inputfiles)
        self._database._u = self
        self.families = {"Base": self._database.base}

        _generate_from_database(self._database.base)

    @property
    def silent(self):
        return self._s

    @silent.setter
    def silent(self, is_silent):
        if not isinstance(is_silent, bool):
            raise TypeError("The value of silence should be True or False.")
        self._s = is_silent


class Time(object):
    @classmethod
    def _mix(cls, other):
        """Creates a time-dependent version of :class:`StructureAttr`."""
        newcls = type(other.__name__, (cls, other), {})
        newcls._derived_class = newcls
        return newcls
