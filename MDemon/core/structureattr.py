import itertools
from abc import ABCMeta, abstractmethod
from bisect import bisect_left

import networkx as nx
import numpy as np

from .. import _ELEMENTS, _STRUCTURE_ATTRS, _STRUCTURE_NAMES, _STRUCTURES
from ..lib.util import asiterable, astuple, flat, iterable, wishnotiterable
from .source import Source1D, Source2D
from .structure import Atom, Bond, Particle, Structure, Topology


def sum_output(func):
    def wrapper(*args, **kwargs):
        original_output = func(*args, **kwargs)
        modified_output = sum(flat(original_output))

        return modified_output

    return wrapper


class SAttrMeta(ABCMeta):
    def __init__(cls, name, bases, classdict):
        ABCMeta.__init__(ABCMeta, name, bases, classdict)
        _STRUCTURE_ATTRS.append(cls)


class StructureAttr(metaclass=SAttrMeta):
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

        database.register_source(self)  # register source and deep source
        database.add_Attr(self)
        self._database = database  # database
        self._timedependent = kwargs.get("timedependent", False)  # time-dependent
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
    def timedependent(self):
        return self._timedependent

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

    @classmethod
    def auto_instancing(
        cls, *valueslist, sid=None, database=None, process_attr=True, **kwargs
    ):
        fname = kwargs.pop("fname", "Base")
        if cls.name in database.attrnames:
            attr = database.__getattribute__(cls.name)
            attr._update_source(valueslist[0], sid=sid, renew=True, **kwargs)
        else:
            attr = cls(*valueslist, sid=sid, database=database, **kwargs)
        if process_attr:
            database.families[fname]._process_attr(attr)
        return attr


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
        if not iterable(numlist[0]):
            valueslist = [np.arange(num) for num in numlist]
        super().__init__(*valueslist, sid=sid, database=database)

    def __getitem__(self, s, sid):
        return wishnotiterable(s._ix)

    @classmethod
    def auto_instancing(
        cls, *valueslist, sid=None, database=None, process_attr=True, **kwargs
    ):
        super().auto_instancing(
            *valueslist, sid=sid, database=database, process_attr=process_attr, **kwargs
        )
        if process_attr:
            fname = sid[0].split("_")[-1]
            database.families[fname].instancing()


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


# TODO: Repair the bug of u1.molecules.mass
class Mass(ParticleAttr):
    name = "mass"
    _dtype = "float"

    @sum_output
    def __getitem__(self, s, sid):
        return super().__getitem__(s, sid)


class Charge(ParticleAttr):
    name = "charge"
    _dtype = "float"


class Temperature(ParticleAttr):
    name = "temperature"
    _dtype = "float"

    @classmethod
    def start_from_velocities(cls, velocities, masses, sid=None, database=None):
        # 常量
        amu_to_kg = 1.66053906660e-27  # 1 amu = 1.66053906660e-27 kg
        angstrom_per_fs_to_m_per_s = 1e5  # 1 Å/fs = 1e5 m/s
        kB = 1.380649e-23  # 玻尔兹曼常数，单位：J/K

        # 单位转换
        velocities_m_per_s = velocities * angstrom_per_fs_to_m_per_s  # 转换速度到 m/s
        masses_kg = masses * amu_to_kg  # 转换质量到 kg

        # 计算每个原子的动能
        kinetic_energies = (
            0.5 * masses_kg[:, np.newaxis] * velocities_m_per_s**2
        )  # 单位：J
        kinetic_energy_per_atom = np.sum(
            kinetic_energies, axis=1
        )  # 对于每个原子，求和 x, y, z 方向的动能

        # 计算温度
        temperatures = (2 * kinetic_energy_per_atom) / (3 * kB)

        return cls(temperatures, sid=sid, database=database)


class AtomAttr(ParticleAttr):
    _sid0 = "Atom_Base"
    name = "atomattr"
    _tclasses = (Atom,)


class Element(AtomAttr):
    name = "element"
    _dtype = "int"

    @classmethod
    def start_from_masses(cls, masses, sid=None, database=None):
        values = np.full((len(masses),), "", dtype=str)
        sorted_elements = _ELEMENTS
        Emasses = [e[0] for e in sorted_elements]

        # 定义一个函数，通过二分搜索找到最接近的元素
        def find_closest_element(mass):
            # 提取所有元素的质量
            pos = bisect_left(Emasses, mass)  # 找到插入点
            if pos == 0:
                return sorted_elements[0][1]
            if pos == len(sorted_elements):
                return sorted_elements[-1][1]
            before = sorted_elements[pos - 1]
            after = sorted_elements[pos]
            # 返回最接近的元素
            if abs(mass - before[0]) < abs(mass - after[0]):
                return before[1]
            else:
                return after[1]

        for i, mass in enumerate(masses):
            values[i] = find_closest_element(mass).number

        cls(values, sid=sid, database=database)


class Valence(AtomAttr):
    name = "valence"
    _dtype = "float"

    @classmethod
    def start_from_bondorders(cls, bondorders, fname="Base", database=None):
        # compomtrx is a csr_matrix
        compomtrx = database.composition._source_register[
            ("Atom_" + fname, "Bond_" + fname)
        ].values

        # 获取CSR矩阵的索引和指针
        indptr = compomtrx.indptr
        indices = compomtrx.indices
        # 使用高级索引和NumPy的add.reduceat来进行高效求和
        values = np.add.reduceat(bondorders[indices], indptr[:-1])

        cls(values, sid=("Atom_" + fname), database=database)


class LonePair(AtomAttr):
    name = "lonepair"
    _dtype = "float"


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

    def _update_source(self, values, sid, N=None, M=None, renew=False, reverse=False):
        # the format of values should be [row,col,data]
        if reverse:
            trsp, sid, values = self._transpose(sid, values)
            N, M = M, N
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

    def union(self, sname0, sname1, process_attr=True):
        """
        atoms and mols are just metaphors.
        """
        G = nx.Graph()
        fname = sname0.split("_")[-1]
        pairs = self.create_pairs(sid=(sname0, sname0))
        G.add_edges_from(pairs)
        components = list(nx.connected_components(G))

        n_mols = len(components)
        n_atoms = self._source_register[(sname0, sname0)].N

        # All about molecule should be updated.
        ixs_mol = np.arange(n_mols, dtype=np.int32)
        Index.auto_instancing(
            ixs_mol,
            sid=(sname1,),
            database=self._database,
            fname=fname,
            process_attr=process_attr,
        )
        # update composition
        row = np.zeros(n_atoms, dtype=np.int32)
        col = np.zeros(n_atoms, dtype=np.int32)
        data = np.ones(n_atoms, dtype=np.float32)
        m = 0
        for i, atmsMol in enumerate(components):
            n_atmsMol = len(atmsMol)
            row[m : m + n_atmsMol] = np.full(n_atmsMol, i, dtype=np.int32)
            col[m : m + n_atmsMol] = list(atmsMol)
            m += n_atmsMol

        compo = Composition.auto_instancing(
            np.array([row, col, data]),
            sid=(sname1, sname0),
            database=self._database,
            N=n_mols,
            M=n_atoms,
            fname=fname,
            process_attr=process_attr,
        )
        if sname0.split("_")[0] != "Atom":
            compo.big_middle_small(
                sname0=sname1, sname1=sname0, sname2="Atom_" + sname0.split("_")[-1]
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
            self._attrnamedic[sid] = cls.abbreviation + "_top"

    def to_connection(self, sname, process_attr=True):
        """
        The default definition of connection is
        sharing the same bond.
        """
        sid = ("Bond_" + sname.split("_")[-1], sname)
        v = self._source_register[sid].values

        # 获取CSR矩阵的行索引和列索引
        rows, cols = v.nonzero()

        # 使用一个字典来收集每一行的非零列索引
        row_dict = {}
        for r, c in zip(rows, cols):
            if r not in row_dict:
                row_dict[r] = []
            row_dict[r].append(c)

        # 收集所有可能的连接对
        connections = []
        for indices in row_dict.values():
            if len(indices) > 1:
                connections.extend(list(itertools.permutations(indices, 2)))

        connections = np.array(list(set(connections)))
        row = connections[:, 0]
        col = connections[:, 1]
        data = np.ones(len(connections), dtype=np.int32)
        Connection.auto_instancing(
            np.array([row, col, data]),
            sid=(sname, sname),
            database=self._database,
            N=v.shape[1],
            M=v.shape[1],
            fname=sname.split("_")[-1],
            process_attr=process_attr,
        )

    def _update_source(
        self, values, sid, N=None, M=None, renew=False, reverse=False, simplemode=False
    ):
        if simplemode:
            # values = [[],[],...,[]]
            n_temp = len(flat(values))
            row = np.zeros(n_temp, dtype=np.int32)
            col = np.zeros(n_temp, dtype=np.int32)
            data = np.ones(n_temp, dtype=np.float32)
            k = 0
            for i, ixs in enumerate(values):
                n_ix = len(ixs)
                row[k : k + n_ix] = np.full(n_ix, i, dtype=np.int32)
                col[k : k + n_ix] = ixs
                k += n_ix
            values = np.array([row, col, data])

        return super()._update_source(values, sid, N, M, renew, reverse)

    def big_middle_small(self, sname0, sname1, sname2):
        v0 = self._source_register[(sname1, sname0)].values
        v2 = self._source_register[(sname1, sname2)].values

        row0, col0 = v0.nonzero()
        middle2big = dict(zip(row0, col0))

        row2, col2 = v2.nonzero()
        n_smalls = len(row2)

        row = np.zeros(n_smalls, dtype=np.int32)
        col = np.zeros(n_smalls, dtype=np.int32)
        dat = np.ones(n_smalls, dtype=np.int32)

        for i, ix in enumerate(row2):
            if ix in middle2big:
                row[i] = middle2big[ix]
                col[i] = col2[i]

        self.auto_instancing(
            np.array([row, col, dat]),
            sid=(sname0, sname2),
            N=v0.shape[1],
            M=v2.shape[1],
            fname=sname0.split("_")[-1],
            database=self._database,
            process_attr=True,
        )


__all__ = []
for attr in _STRUCTURE_ATTRS:
    __all__.append(attr.__name__)
