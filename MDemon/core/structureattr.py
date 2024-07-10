import itertools
from abc import ABCMeta, abstractmethod

import networkx as nx
import numpy as np
from scipy import sparse

from .. import _STRUCTURE_ATTRS, _STRUCTURE_NAMES, _STRUCTURES
from ..lib.util import asiterable, astuple, wishnotiterable
from .source import Source1D, Source2D
from .structure import Atom, Bond, Molecule, Particle, Structure, Topology


class SAttrMeta(ABCMeta):
    def __init__(cls, name, bases, classdict):
        ABCMeta.__init__(ABCMeta, name, bases, classdict)
        _STRUCTURE_ATTRS.append(cls)


class StructureAttr(object, metaclass=SAttrMeta):
    """
    Base class of structure attributes.

    Parameters
    ----------
    valueslist : `list`
        Initial input values when creating a new \
        :class:`StructureAttr` object, which will \
        be converted to the first `dict.value` of \
        `source` after initialization.
    sid : `tuple`
        The structural id of input values \
        as well as the first `dict.key` in `source`.

    Attributes
    ----------
    source : `dict`
        Core data of :class:`StructureAttr` object.
    timedependent : `bool`
        Whether the :class:`StructureAttr` is time-dependent \
        or not. If `timedependent` is `False`, \
        the `source` will not change with the global timestep.
    timestep : `int`
        The timestep inside the attribute. When you \
        change the global timestep, `ts` will \
        not change automatically until using \
        :function:`__getitem__` to access `source`.
    deepsource : `list`
        The data depository for time-dependent \
        :class:`StructureAttr` objects.
    database : `Database`
        The pointer to its database.
    name : `str`
        The name of :class:`StructureAttr`.
    _sid0 : `tuple`
        The default structural id of input values.
    _dtype : `str`
        Should be 'bool', 'int' or 'float'

    """

    name = "structureattr"  # name
    _sid0 = ("Structure_Base",)  # default structural base
    _tclasses = (Structure,)  # targeted classes
    _dtype = ""  # data type

    def __init__(self, *valueslist, sid=None, database=None, **kwargs):
        """
        We register only one structural id during the \
        initialization considering of the readability of codes. \
        When :class:`StructureAttr` is time-dependent, \
        a `deepsource` will be created for convenient data \
        conversions between different timesteps. \
        While in the time-independent situation, \
        we just skip the establishment of `deepsource`.

        Parameters
        ----------
        valueslist : `list`
            Initial input values when creating a new \
            :class:`StructureAttr` object, which will \
            be converted to the first `dict.value` of \
            `source` after initialization.
        sid : `tuple`
            The structural id of input values \
            as well as the first `dict.key` in `source`.
        """

        database.register_source(self)  # source
        database.register_deep_source(self)  # deep source
        database.add_Attr(self)
        self._database = database  # database
        self.timedependent = False  # time-dependent
        self.timestep = 0  # timestep

        if sid is None:
            sid = self._sid0
        else:
            sid = astuple(sid)

        if self.timedependent:
            self._update_deepsource(valueslist, sid)
            self.source = self.deepsource[0]
        else:
            self._update_source(valueslist[0], sid, **kwargs)

    @property
    def targeted_classes(self):
        tc = []
        for cls in _STRUCTURES:
            if issubclass(cls, self._tclasses):
                tc.append(cls)
        return tc

    def _update_deepsource(self, valueslist, sid):
        """
        Add a valueslist to `deepsource`.

        Parameters
        ----------
        valueslist : `list`
            Input values in the format of `list`.
        sid : `tuple`
            The structural id of input values.
        """
        for i, values in enumerate(valueslist):
            values = self._update_source(values, sid)
            try:
                self.deepsource[i][sid] = values
            except IndexError:
                self.deepsource.append({sid: values})

    @abstractmethod
    def _update_source(self, **kwargs):
        """
        Convert the format of input values into the pre-set format. \
        If `ix` is not `None`, update `source`. \
        If `ix` equal -1, update the whole `source`.

        Parameters
        ----------
        values : `vector` or `matrix`
            The order of value should be consistent \
            with the order of ix.
        sid : `tuple`
            The structural id of input values.
        ix : `numpy.ndarray` or `numpy.intp`
            The index of values.
        Returns
        -------
        values : `numpy.ndarray` or `scipy.sparse.csr_matrix`
        """

    def __getitem__(self, s, sid):
        """
        We always call `__getitem__` in the way of \
        the `Structure` object's attributes, which \
        means, we wrap the `StructureAttr` object \
        into a `property` attribute of targeting \
        `Structure` objects. If no corresponding \
        sid is found in `self.source`, we will \
        try to access those values in a "basic type version".

        Parameters
        ----------
        s : `Structure`
            The `Structure` object who calls this function.
        sid : `tuple`
            The structural id of the values you want to access.
        Returns
        -------
        values : `numpy.ndarray` or `scipy.sparse.csr_matrix`
        """
        if self.timedependent:
            self._reload_source()
        try:
            values = self._parse_source(sid, s.ix, s.detailed)
            return wishnotiterable(values)
        except KeyError:
            sid, ixs = self._basic_type_version(s, sid)
            return np.asarray([self._parse_source(sid, ix, s.detailed) for ix in ixs])

    def _reload_source(self):
        """
        Reload `source` if the global timestep changed.
        """
        u = self._database._u
        dt = u.timestep - self.timestep
        if dt != 0:
            sign = np.sign(dt)
            t1 = self.timestep + sign
            t2 = dt + self.timestep + sign
            for i in range(t1, t2, sign):
                vdic = self.deepsource[i]
                self.source += vdic

    @abstractmethod
    def _parse_source(self, **kwargs):
        """
        We use this function to transform \
        the storage format to a more \
        readable format.

        Parameters
        ----------
        values : `numpy.ndarray` or `scipy.sparse.csr_matrix`
        ix : `numpy.ndarray` or `numpy.intp`
        detailed : `bool`

        Returns
        -------
        values : `numpy.ndarray` or `dict`
        """

    def _basic_type_version(self, s, sid):
        """
        Convert the original sid into a basic type version.

        Parameters
        ----------
        s : `Structure`
            The `Structure` object who calls this function.
        sid : `tuple`
            The structural id of the values you want to access.
        Returns
        -------
        sid : `tuple`
            The basic type version of original sid.
        ixs : `numpy.ndarray`
            The indices of corresponding basic objects.
        """
        sid = list(sid)
        if isinstance(s, Particle):
            ixs = s.atms
            sid[0] = "Atom_" + s._fname
        elif isinstance(s, Topology):
            ixs = s.bnds
            sid[0] = "Bond_" + s._fname
        sid = tuple(sid)
        return sid, ixs

    def __setitem__(self, s, sid, values):
        """
        In current version, we always set the values \
        in a holistic way, which means we do not offer \
        the option to select indices manually.

        Parameters
        ----------
        s : `Structure`
            The `Structure` object who calls this function.
        sid : `tuple`
            The structural id of the values you want to access.
        values : `vector` or `matrix`
            You should provide values for all indices in `s`.
        """
        if sid in self.source:
            self._update_source(values, sid, s.ix)

    @abstractmethod
    def import_attrname(self, **kwargs):
        """
        When we attach a `StructureAttr` object to \
        a subclass of :class:`Structure`, we need to \
        tell the subclass which `attrname` it need \
        have to access corresponding values.

        Parameters
        ----------
        sid : `tuple`

        Returns
        -------
        attrname : `str`
        """

    @abstractmethod
    def import_sids(self, **kwargs):
        """
        Sometimes there are more than one `sid` \
        for certain :class:`Structure` subclass \
        in an `StructureAttr` object. As a result, \
        we need import all related `sids`.

        Parameters
        ----------
        cls : `Structure` class

        Returns
        -------
        sids : `tuple` stored in `list`
        """

    @property
    def snames(self):
        """
        Return all `snames` mentioned in all `sids`.

        Returns
        -------
        snames : `tuple`
        """
        snames = []
        for sid in self.source:
            snames.extend(list(sid))
        return tuple(set(snames))

    @classmethod
    def _subclass(cls, **kwargs):
        newcls = type(cls.__name__, (cls,), kwargs)

        return newcls


class StructureAttr1D(StructureAttr):
    """
    :class:`StructureAttr1D` contains all structure \
    attributes with len(`sid`) = 1.
    """

    name = "structureattr1D"

    def _update_source(self, values, sid, ix=None, renew=False):
        values = np.asarray(values)
        if sid in self._source_register and not renew:
            source = self._source_register[sid]
            valix = [ix, values]
            source.values = valix
        else:
            source = Source1D(self._dtype, values)
            self._source_register[sid] = source

    def _parse_source(self, sid, ix, detailed):
        source = self._source_register[sid]
        return source.values[ix]

    def import_attrname(self, *args):
        return self.name

    @staticmethod
    def import_sids(cls):
        if isinstance(cls, type):
            return [(cls.import_sname(),)]
        else:
            return [(cls.sname,)]


class Absence(StructureAttr1D):
    name = "absence"
    _dtype = "bool"


class Silence(StructureAttr1D):
    name = "silence"
    _dtype = "bool"


class Freeze(StructureAttr1D):
    name = "freeze"
    _dtype = "bool"


class ID(StructureAttr1D):
    name = "id"
    _dtype = "int"


class Index(StructureAttr1D):
    name = "ix"
    _dtype = "int"

    def __init__(self, *numlist, sid, database):
        valueslist = [np.arange(num) for num in numlist]
        super().__init__(*valueslist, sid=sid, database=database)

    def __getitem__(self, s, sid):
        return wishnotiterable(s._ix)


class Species(StructureAttr1D):
    name = "species"
    _dtype = "int"


class Coordinate(StructureAttr1D):
    name = "coordinate"
    _dtype = "float"


class Velocity(StructureAttr1D):
    name = "velocity"
    _dtype = "float"


class ParticleAttr(StructureAttr1D):
    _sid0 = "Particle_Base"
    name = "particleattr"
    _tclasses = (Particle,)


class Mass(ParticleAttr):
    name = "mass"
    _dtype = "float"


class Charge(ParticleAttr):
    name = "charge"
    _dtype = "float"


class AtomAttr(ParticleAttr):
    _sid0 = "Atom_Base"
    name = "atomattr"
    _tclasses = (Atom,)


class Element(AtomAttr):
    name = "element"
    _dtype = str


class TopologyAttr(StructureAttr1D):
    _sid0 = "Topology_Base"
    name = "topologyattr"
    _tclasses = (Topology,)


class BondAttr(TopologyAttr):
    _sid0 = "Bond_Base"
    name = "bondattr"
    _tclasses = (Bond,)


class BondOrder(BondAttr):
    name = "bondorder"
    _dtype = "float"


class StructureAttr2D(StructureAttr, metaclass=ABCMeta):
    """
    :class:`StructureAttr1D` contains all structure \
    attributes with len(`sid`) = 2.
    """

    _attrnamedic = {}

    def import_attrname(self, sid):
        return self._attrnamedic[sid]

    def import_sids(self, cls):
        sids = []
        if isinstance(cls, type):
            sname = cls.import_sname()
        else:
            sname = cls.sname

        for sid in self._source_register:
            if sname == sid[0]:
                sids.append(sid)
        return sids

    def _parse_source(self, sid, ix, detailed):
        v = self._source_register[sid].values
        ix = asiterable(ix)
        values = []
        for i in ix:
            indices = v[i].indices
            if detailed:
                data = v[i].data
                values.append(dict(zip(indices, data)))
            else:
                values.append(indices)
        return wishnotiterable(values)

    def _update_source(self, values, sid, N=None, M=None, renew=False):
        # the format of values should be [row,col,data]
        values = np.array(values)
        if renew:
            self._init_mtrx(sid, values, N, M)
        else:
            try:
                self._add(sid, values)
            except KeyError:
                self._init_mtrx(sid, values, N, M)

    def _init_mtrx(self, sid, values, N, M):
        """
        Initiate the value matrix, then register \
        `attrname` for every single `sid`.

        Parameters
        ----------
        sid : `tuple`
        values : `scipy.sparse.csr_matrix`
        """
        source = Source2D(self._dtype, values, N=N, M=M)
        self._source_register[sid] = source
        self._register_attrname(sid, 0)
        trsp, sid1, values1 = self._transpose(sid, values)
        if trsp:
            source1 = Source2D(self._dtype, values1, N=M, M=N)
            self._source_register[sid1] = source1
            self._register_attrname(sid1, 1)

    @abstractmethod
    def _register_attrname(self, **kwargs):
        """
        Different from 1D attr, the attrname in \
        2D attr changes with sid.

        Parameters
        ----------
        sid : `tuple`
        """

    def _add(self, sid, values):
        source = self._source_register[sid]
        source.values = values
        trsp, sid1, values1 = self._transpose(sid, values)
        if trsp:
            source1 = self._source_register[sid1]
            source1.values = values1

    @staticmethod
    def _transpose(sid, values):
        sid1 = (sid[1], sid[0])
        trsp = sid1 != sid
        if trsp:
            values = np.array([values[1], values[0], values[2]])
        return trsp, sid1, values

    def create_pairs(self, **kwargs):
        sid = kwargs.get("sid", ("Atom_Base", "Atom_Base"))
        v = self._source_register[sid].values

        if sid[0] == sid[1]:
            selfcombi = True
        else:
            selfcombi = False

        row_indices, col_indices = v.nonzero()

        if selfcombi:
            mask = row_indices <= col_indices
            unique_row_indices = row_indices[mask]
            unique_col_indices = col_indices[mask]
            unique_indices = np.vstack((unique_row_indices, unique_col_indices)).T
            return unique_indices

        return np.vstack((row_indices, col_indices)).T


class Connection(StructureAttr2D):
    name = "connection"
    _dtype = "int"

    def _register_attrname(self, sid, x):
        self._attrnamedic[sid] = "neighbors"

    def to_molecules(self, fname="Base"):
        G = nx.Graph()
        sname = "Atom_" + fname
        pairs = self.create_pairs(sid=(sname, sname))
        G.add_edges_from(pairs)
        components = list(nx.connected_components(G))
        n_mols = len(components)
        n_atoms = self._source_register[(sname, sname)].N

        # All about molecule should be updated.
        mol_ids = np.arange(1, n_mols + 1, dtype=np.int32)
        mol_ixs = np.arange(n_mols, dtype=np.int32)
        msname = "Molecule_" + fname
        self._database.id._update_source(mol_ids, (msname,), renew=True)
        self._database.ix._update_source(mol_ixs, (msname,), renew=True)

        # uodate composition
        row = np.zeros(n_atoms, dtype=np.int32)
        col = np.zeros(n_atoms, dtype=np.int32)
        data = np.zeros(n_atoms, dtype=np.float32)
        m = 0
        for i, atmsMol in enumerate(components):
            n_atmsMol = len(atmsMol)
            row[m : m + n_atmsMol] = np.full(n_atmsMol, i, dtype=np.int32)
            col[m : m + n_atmsMol] = list(atmsMol)
            m += n_atmsMol
        self._database.composition._update_source(
            np.array([row, col, data]), (msname, sname), N=n_mols, M=n_atoms, renew=True
        )


class Composition(StructureAttr2D):
    """
    The squence of two snames in sid is very important. \
    We preset that the second sname is the subset of the \
    first sname, for example: ('Molecule_Base','Atom_Base') \
    is a correct input sid while ('Atom_Base','Molecule_Base') \
    is a wrong one.
    """

    name = "composition"
    _dtype = "float"

    def _register_attrname(self, sid, x):
        cls = _STRUCTURE_NAMES[sid[1].split("_")[0]]
        if x == 0:
            self._attrnamedic[sid] = cls.abbreviation + "s"
        elif x == 1:
            self._attrnamedic[sid] = cls.abbreviation + "_ix"

    def to_connection(self, sname):
        sid = (
            sname,
            "Atom_" + sname.split("_")[-1],
        )
        v = self._source_register[sid].values

        row = np.zeros(0, dtype=np.int32)
        col = np.zeros(0, dtype=np.int32)
        data = np.zeros(0, dtype=np.int32)
        for i in range(v.shape[0]):
            list_ = v[i].indices
            if len(list_) > 1:
                combinations = np.array(list(itertools.permutations(list_, 2)))
                row_ = combinations[:, 0]
                col_ = combinations[:, 1]
                data_ = np.full(len(row_), i, dtype=np.int32)
                row = np.concatenate((row, row_))
                col = np.concatenate((col, col_))
                data = np.concatenate((data, data_))
        return np.array([row, col, data])


__all__ = []
for attr in _STRUCTURE_ATTRS:
    __all__.append(attr.__name__)
