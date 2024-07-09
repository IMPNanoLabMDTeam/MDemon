import numpy as np

from .. import _STRUCTURES
from .structureattr import Absence, Freeze, Index, Silence


class Database(object):
    _u = None

    def __init__(self, n_atoms, **kwargs):
        self._source_register = {}
        self._deep_source_register = {}
        self.attrs = []

        self._n_atoms = n_atoms

        self.time_dependent = kwargs.get("time_dependent", False)
        self.set_absence = kwargs.get("set_absence", False)
        self.set_silence = kwargs.get("set_silence", False)
        self.set_freeze = kwargs.get("set_freeze", False)

        Index(n_atoms, sid=("Atom_Base",), database=self)
        self._base = Family(database=self)

    @property
    def time_dependent(self):
        return self._time_dependent

    @time_dependent.setter
    def time_dependent(self, val):
        if val and not self._time_dependent:
            self._time_dependent = True
            self._deep_source_register = {}
        else:
            self._time_dependent = False

    @property
    def set_absence(self):
        return self._set_absence

    @set_absence.setter
    def set_absence(self, val):
        if val and not self._set_absence:
            self._set_absence = True
            attr = Absence(
                np.array([False for _ in range(self.n_atoms)]),
                sid=("Atom_Base",),
                database=self,
            )
            self.add_Attr(attr)
        else:
            self._set_absence = False

    @property
    def set_silence(self):
        return self._set_silence

    @set_silence.setter
    def set_silence(self, val):
        if val and not self._set_silence:
            self._set_silence = True
            attr = Silence(
                np.array([False for _ in range(self.n_atoms)]),
                sid=("Atom_Base",),
                database=self,
            )
            self.add_Attr(attr)
        else:
            self._set_silence = False

    @property
    def set_freeze(self):
        return self._set_freeze

    @set_freeze.setter
    def set_freeze(self, val):
        if val and not self._set_freeze:
            self._set_freeze = True
            attr = Freeze(
                np.array([False for _ in range(self.n_atoms)]),
                sid=("Atom_Base",),
                database=self,
            )
            self.add_Attr(attr)
        else:
            self._set_freeze = False

    def add_Attr(self, attr):
        """Add a new Attr to the Database.

        Parameters
        ----------
        attr : :class:`StructureAttr` or :class:`UniverseAttr`

        """
        self.attrs.append(attr)
        attr._database = self
        self.__setattr__(attr.name, attr)

    @property
    def n_atoms(self):
        return self._n_atoms

    @property
    def base(self):
        return self._base

    def register_source(self, attr):
        self._source_register[attr.name] = {}
        attr._source_register = self._source_register[attr.name]

    def register_deep_source(self, attr):
        if self.time_dependent:
            self._deep_source_register[attr.name] = {}
            attr._deep_source_register = self._deep_source_register[attr.name]


class Family(object):
    def __init__(self, fname="Base", database=None):
        self.name = fname
        self._database = database

    def instancing(self):
        for cls in _STRUCTURES:
            cls1 = self._u.families[self.name]._classes[cls]
            sid = (cls1.import_sname(),)
            attrname = cls1.__name__.lower() + "s"
            index = self._database.ix
            if sid in index._source_register:
                n = len(index._source_register[sid].values)
                s = cls1(np.arange(n), self._u)
                self.__setattr__(attrname, s)
                if self.name == "Base":
                    self._u.__setattr__(attrname, s)

    def _process_attr(self, attr):
        """Squeeze a structureattr for its information

        Grabs:
         - Structure properties (attribute access)
         - Transplant methods (unaccessible now)
        """

        for cls in attr.targeted_classes:
            self._classes[cls]._add_prop(attr)

    @property
    def _u(self):
        return self._database._u
