from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import sparse

from .. import _STRUCTURE_NAMES, _STRUCTURES
from ..lib.util import asiterable, astuple, wishnotiterable
from .structure import Atom, Bond, Molecule, Particle, Structure, Topology


class StructureAttr(object, metaclass=ABCMeta):
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
    _dtype : `int` or `float`
        The data type of `source`.

    """

    name = "structureattr"  # name
    _sid0 = ("Structure_Base",)  # default structural base
    _tclasses = (Structure,)  # targeted classes
    _dtype = None  # data type

    def __init__(self, *valueslist, sid=None):
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
        self.source = {}  # source
        self.deepsource = []  # deep source
        self._database = None  # database
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
            self._update_source(valueslist[0], sid, -1)

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
            ixs = s.atm_ix
            sid[0] = s.family.atom.sname
        elif isinstance(s, Topology):
            ixs = s.bnd_ix
            sid[0] = s.family.bond.sname
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


class StructureAttr1D(StructureAttr):
    """
    :class:`StructureAttr1D` contains all structure \
    attributes with len(`sid`) = 1.
    """

    name = "structureattr1D"

    def _update_source(self, values, sid, ix=None):
        values = np.asarray(values, dtype=self._dtype)
        if ix:
            if ix != -1:
                self.source[sid][ix] = values
            else:
                self.source[sid] = values
        return {sid: values}

    def _parse_source(self, sid, ix, detailed):
        values = self.source[sid]
        return values[ix]

    def import_attrname(self, *args):
        return self.name

    def import_sids(self, cls):
        if isinstance(cls, type):
            return [(cls.import_sname(),)]
        else:
            return [(cls.sname,)]


class Existence(StructureAttr1D):
    name = "existence"
    _dtype = bool


class Silence(StructureAttr1D):
    name = "silence"
    _dtype = bool


class Freeze(StructureAttr1D):
    name = "freeze"
    _dtype = bool


class ID(StructureAttr1D):
    name = "id"
    _dtype = np.intp


class Index(StructureAttr1D):
    name = "ix"
    _dtype = np.intp

    def __init__(self, *numlist, sid):
        valueslist = [np.arange(num) for num in numlist]
        super().__init__(valueslist, sid=sid)

    def __getitem__(self, s, sid):
        return wishnotiterable(s._ix)


class Species(StructureAttr1D):
    name = "species"
    _dtype = np.intp


class Coordinate(StructureAttr1D):
    name = "coordinate"
    _dtype = np.float32


class Velocity(StructureAttr1D):
    name = "velocity"
    _dtype = np.float32


class ParticleAttr(StructureAttr1D):
    _sid0 = "Particle_Base"
    name = "particleattr"
    _tclasses = (Particle,)


class Mass(ParticleAttr):
    name = "mass"
    _dtype = np.float32


class Charge(ParticleAttr):
    name = "charge"
    _dtype = np.float32


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
    _dtype = np.float32


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

        for sid in self.source:
            if sname == sid[0]:
                sids.append(sid)
        return sids

    def _parse_source(self, sid, ix, detailed):
        s = self.source[sid]
        ix = asiterable(ix)
        values = []
        for i in ix:
            indices = s[i].indices
            if detailed:
                data = s[i].data
                values.append(dict(zip(indices, data)))
            else:
                values.append(indices)
        return values

    # TODO: Manually modify values with ix.
    def _update_source(self, values, sid, ix):
        values = sparse.csr_matrix(values, dtype=self._dtype)
        if ix is not None:
            if sid not in self.source:
                self._init_mtrx(sid, values)
            else:
                self._add(sid, values)

    def _init_mtrx(self, sid, vlmtrx):
        """
        Initiate the value matrix, then register \
        `attrname` for every single `sid`.

        Parameters
        ----------
        sid : `tuple`
        vlmtrx : `scipy.sparse.csr_matrix`
        """
        self.source[sid] = sparse.csr_matrix(vlmtrx.shape, dtype=self._dtype)
        self._register_attrname(sid, 0)
        self.source[sid] += vlmtrx
        trsp, sid1, vlmtrx1 = self._transpose(sid, vlmtrx)
        if trsp:
            self.source[sid1] = sparse.csr_matrix(vlmtrx1.shape, dtype=self._dtype)
            self._register_attrname(sid1, 1)
            self.source[sid1] += vlmtrx1

    @abstractmethod
    def _register_attrname(self, **kwargs):
        """
        Different from 1D attr, the attrname in \
        2D attr changes with sid.

        Parameters
        ----------
        sid : `tuple`
        """

    def _add(self, sid, vlmtrx):
        self.source[sid] += vlmtrx
        trsp, sid1, vlmtrx1 = self._transpose(sid, vlmtrx)
        if trsp:
            self.source[sid1] += vlmtrx1

    @staticmethod
    def _transpose(sid, vlmtrx):
        sid1 = (sid[1], sid[0])
        trsp = sid1 != sid
        if trsp:
            vlmtrx = np.transpose(vlmtrx)
        return trsp, sid1, vlmtrx


class Connection(StructureAttr2D):
    def _register_attrname(self, sid, x):
        self._attrnamedic[sid] = "neighbors"


class Composition(StructureAttr2D):
    """
    The squence of two snames in sid is very important. \
    We preset that the second sname is the subset of the \
    first sname, for example: ('Molecule_Base','Atom_Base') \
    is a correct input sid while ('Atom_Base','Molecule_Base') \
    is a wrong one.
    """

    def _register_attrname(self, sid, x):
        cls = _STRUCTURE_NAMES[sid[1].split("_")[0]]
        if x == 0:
            self._attrnamedic[sid] = cls.abbreviation + "s"
        elif x == 1:
            self._attrnamedic[sid] = cls.abbreviation + "_ix"
