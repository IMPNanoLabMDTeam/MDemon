import numpy as np

from .. import _STRUCTURES
from .structureattr import Existence, Freeze, Index, Silence


class Database(object):
    _u = None

    def __init__(self, ids, attrs: list):
        atm = "Atom_Base"
        aid = (atm,)
        n_atoms = len(ids.source[aid])
        ixs = Index(n_atoms, sid=aid)
        for cls in _STRUCTURES:
            sid = (cls.import_sname(),)
            if sid in ids.source:
                n = len(ids.source[sid])
                ixs._update_source(n, sid=sid)

        exists = Existence(np.array([True for _ in range(n_atoms)]), sid=aid)
        silens = Silence(np.array([False for _ in range(n_atoms)]), sid=aid)
        frees = Freeze(np.array([False for _ in range(n_atoms)]), sid=aid)

        attrs.extend((ixs, exists, silens, frees))
        # attach the StructureAttrs
        self.attrs = []
        for attr in attrs:
            self.add_StructureAttr(attr)

        self.base = Family(database=self)

    def add_StructureAttr(self, attr):
        """Add a new StructureAttr to the Database.

        Parameters
        ----------
        attr : :class:`StructureAttr`

        """
        self.attrs.append(attr)
        attr._database = self
        self.__setattr__(attr.name, attr)


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
            if sid in index.source:
                n = len(index.source[sid])
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
