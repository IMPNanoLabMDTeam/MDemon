import numpy as np

from ..core.database import Database
from ..core.structureattr import (
    ID,
    Charge,
    Composition,
    Coordinate,
    Element,
    Mass,
    Species,
    Temperature,
    Velocity,
)
from ..core.universe import Box
from .base import ReaderBase, squash_by

# Sections will all start with one of these words
# and run until the next section title
SECTIONS = set(
    [
        "Atoms",  # Molecular topology sections
        "Velocities",
        "Masses",
        "Ellipsoids",
        "Lines",
        "Triangles",
        "Bodies",
        "Bonds",  # Forcefield sections
        "Angles",
        "Dihedrals",
        "Impropers",
        "Pair",
        "Pair LJCoeffs",
        "PairIJ Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
        "BondBond Coeffs",  # Class 2 FF sections
        "BondAngle Coeffs",
        "MiddleBondTorsion Coeffs",
        "EndBondTorsion Coeffs",
        "AngleTorsion Coeffs",
        "AngleAngleTorsion Coeffs",
        "BondBond13 Coeffs",
        "AngleAngle Coeffs",
    ]
)
# We usually check by splitting around whitespace, so check
# if any SECTION keywords will trip up on this
# and add them
for val in list(SECTIONS):
    if len(val.split()) > 1:
        SECTIONS.add(val.split()[0])


HEADERS = set(
    [
        "atoms",
        "bonds",
        "angles",
        "dihedrals",
        "impropers",
        "atom types",
        "bond types",
        "angle types",
        "dihedral types",
        "improper types",
        "extra bond per atom",
        "extra angle per atom",
        "extra dihedral per atom",
        "extra improper per atom",
        "extra special per atom",
        "ellipsoids",
        "lines",
        "triangles",
        "bodies",
        "xlo xhi",
        "ylo yhi",
        "zlo zhi",
        "xy xz yz",
    ]
)


class DATAReader(ReaderBase):
    """ """

    format = "DATA"

    def iterdata(self):
        with open(self.filename) as f:
            for line in f:
                line = line.partition("#")[0].strip()
                if line:
                    yield line

    def grab_datafile(self):
        """Split a data file into dict of header and sections

        Returns
        -------
        header - dict of header section: value
        sections - dict of section name: content
        """
        f = list(self.iterdata())

        starts = [i for i, line in enumerate(f) if line.split()[0] in SECTIONS]
        starts += [None]

        header = {}
        for line in f[: starts[0]]:
            for token in HEADERS:
                if line.endswith(token):
                    header[token] = line.split(token)[0]
                    continue

        sects = {f[l]: f[l + 1 : starts[i + 1]] for i, l in enumerate(starts[:-1])}

        return header, sects

    @staticmethod
    def _interpret_atom_style(atom_style):
        """Transform a string description of atom style into a dict

        Required fields: id, type, x, y, z
        Optional fields: resid, charge

        eg: "id resid type charge x y z"
        {'id': 0,
         'resid': 1,
         'type': 2,
         'charge': 3,
         'x': 4,
         'y': 5,
         'z': 6,
        }
        """
        style_dict = {}

        atom_style = atom_style.split()

        for attr in ["id", "type", "resid", "charge", "x", "y", "z"]:
            try:
                location = atom_style.index(attr)
            except ValueError:
                pass
            else:
                style_dict[attr] = location

        reqd_attrs = ["id", "type", "x", "y", "z"]
        missing_attrs = [attr for attr in reqd_attrs if attr not in style_dict]
        if missing_attrs:
            raise ValueError(
                "atom_style string missing required field(s): {}"
                "".format(", ".join(missing_attrs))
            )

        return style_dict

    def parse(self, **kwargs):
        """Parses a LAMMPS_ DATA file.

        Returns
        -------
        MDAnalysis Topology object.
        """
        # Can pass atom_style to help parsing
        try:
            self.style_dict = self._interpret_atom_style(kwargs["atom_style"])
        except KeyError:
            self.style_dict = None

        head, sects = self.grab_datafile()

        try:
            masses = self._parse_masses(sects["Masses"])
        except KeyError:
            masses = None

        if "Atoms" not in sects:
            raise ValueError("Data file was missing Atoms section")

        # try:
        ids, species, compo, dbase, order, masses = self._parse_atoms(
            sects["Atoms"], masses
        )
        # except Exception:
        #     errmsg = (
        #         "Failed to parse atoms section.  You can supply a description "
        #         "of the atom_style as a keyword argument, "
        #         "eg mda.Universe(..., atom_style='id resid x y z')"
        #     )
        #     raise ValueError(errmsg) from None

        # create mapping of id to index (ie atom id 10 might be the 0th atom)

        sname_atm = "Atom_Base"
        mapping = {
            atom_id: i
            for i, atom_id in enumerate(ids._source_register[("Atom_Base",)].values)
        }
        n_atoms = len(mapping)
        if "Velocities" in sects:
            velocities = self._parse_vel(sects["Velocities"], order)
        else:
            velocities = np.zeros((n_atoms, 3))

        Velocity(velocities, sid=(sname_atm,), database=dbase)
        Temperature.start_from_velocities(
            velocities, masses, sid=(sname_atm,), database=dbase
        )

        for sname, L, nentries in [
            ("Bond_Base", "Bonds", 2),
            ("Angle_Base", "Angles", 3),
            ("Dihedral_Base", "Dihedrals", 4),
            ("Improper_Base", "Impropers", 4),
        ]:
            try:
                id_, type, sect = self._parse_bond_section(sects[L], nentries, mapping)
                sid = (sname,)
                ids._update_source(id_, sid)
                dbase.ix._update_source(np.arange(len(id_)), (sname,))
                species._update_source(type, sid)
                mtrx, N, M = self.bondsect2mtrx(sect, n_atoms, nentries)

                if sname == "Bond_Base":
                    compo._update_source(mtrx, (sname, sname_atm), N, M, reverse=True)
                    compo.to_connection(sname_atm, process_attr=False)
                else:
                    compo._update_source(mtrx, (sname, sname_atm), N, M)
            except KeyError:
                id_, type, sect = [], [], []
        Box(self._parse_box(head), database=dbase)
        return dbase

    @staticmethod
    def bondsect2mtrx(sect, n_atoms, nentries):
        N = len(sect)
        M = n_atoms
        data = np.array([1 for _ in range(N * nentries)])
        row = np.array([i for i in range(N) for _ in range(nentries)])
        col = []
        for i in sect:
            col.extend(list(i))
        col = np.asarray(col)
        values = np.array([row, col, data])
        return values, N, M

    def read_DATA_timestep(self, n_atoms, TS_class, TS_kwargs, atom_style=None):
        """Read a DATA file and try and extract x, v, box.

        - positions
        - velocities (optional)
        - box information

        Fills this into the Timestep object and returns it

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.18.0
           Added atom_style kwarg
        """
        if atom_style is None:
            self.style_dict = None
        else:
            self.style_dict = self._interpret_atom_style(atom_style)

        header, sects = self.grab_datafile()

        unitcell = self._parse_box(header)

        try:
            positions, ordering = self._parse_pos(sects["Atoms"])
        except KeyError as err:
            errmsg = f"Position information not found: {err}"
            raise OSError(errmsg) from None

        if "Velocities" in sects:
            velocities = self._parse_vel(sects["Velocities"], ordering)
        else:
            velocities = None

        ts = TS_class.from_coordinates(positions, velocities=velocities, **TS_kwargs)
        ts.dimensions = unitcell

        return ts

    def _parse_pos(self, datalines):
        """Strip coordinate info into np array"""
        pos = np.zeros((len(datalines), 3), dtype=np.float32)
        # TODO: could maybe store this from topology parsing?
        # Or try to reach into Universe?
        # but ugly because assumes lots of things, and Reader should be standalone
        ids = np.zeros(len(pos), dtype=np.int32)

        if self.style_dict is None:
            n = len(datalines[0].split())
            if n in (5, 8):
                # atomic format: id type x y z (+ flag1 flag2 flag3)
                style_dict = {"id": 0, "x": 2, "y": 3, "z": 4}
            elif n in (6, 9):
                # molecular format: id resid type x y z (+ flag1 flag2 flag3)
                style_dict = {"id": 0, "x": 3, "y": 4, "z": 5}
            elif n in (7, 10):
                # full format: id resid type charge x y z (+ flag1 flag2 flag3)
                style_dict = {"id": 0, "x": 4, "y": 5, "z": 6}
            else:
                # fallback
                style_dict = {"id": 0, "x": 3, "y": 4, "z": 5}
        else:
            style_dict = self.style_dict

        for i, line in enumerate(datalines):
            line = line.split()

            ids[i] = line[style_dict["id"]]

            pos[i, :] = [
                line[style_dict["x"]],
                line[style_dict["y"]],
                line[style_dict["z"]],
            ]

        order = np.argsort(ids)
        pos = pos[order]

        # return order for velocities
        return pos, order

    def _parse_vel(self, datalines, order):
        """Strip velocity info into np array

        Parameters
        ----------
        datalines : list
          list of strings from file
        order : np.array
          array which rearranges the velocities into correct order
          (from argsort on atom ids)

        Returns
        -------
        velocities : np.ndarray
        """
        vel = np.zeros((len(datalines), 3), dtype=np.float32)

        for i, line in enumerate(datalines):
            line = line.split()
            vel[i] = line[1:4]

        vel = vel[order]

        return vel

    def _parse_bond_section(self, datalines, nentries, mapping):
        """Read lines and strip information

        Arguments
        ---------
        datalines : list
          the raw lines from the data file
        nentries : int
          number of integers per line
        mapping : dict
          converts atom_ids to index within topology

        Returns
        -------
        types : tuple of strings
          type of the bond/angle/dihedral/improper
        indices : tuple of ints
          indices of atoms involved
        """
        section = []
        type = []
        id_ = []
        for line in datalines:
            line = line.split()
            # map to 0 based int
            section.append(tuple([mapping[int(x)] for x in line[2 : 2 + nentries]]))
            type.append(line[1])
            id_.append(line[0])
        return tuple(id_), tuple(type), tuple(section)

    def _parse_atoms(self, datalines, massdict=None):
        """Creates a Topology object

        Adds the following attributes
         - resid
         - type
         - masses (optional)
         - charge (optional)

        Lammps atoms can have lots of different formats,
        and even custom formats

        http://lammps.sandia.gov/doc/atom_style.html

        Treated here are
        - atoms with 7 fields (with charge) "full"
        - atoms with 6 fields (no charge) "molecular"

        Arguments
        ---------
        datalines - the relevent lines from the data file
        massdict - dictionary relating type to mass

        Returns
        -------
        top - Topology object
        """
        n_atoms = len(datalines)

        if self.style_dict is None:
            # Fields per line
            n = len(datalines[0].split())
            if n in (5, 8):
                # atomic format: id type x y z (+ 可能的flag1 flag2 flag3)
                # 5字段: id type x y z
                # 8字段: id type x y z flag1 flag2 flag3
                sd = {"id": 0, "type": 1, "coord": 2}
            elif n in (6, 9):
                # molecular format: id resid type x y z (+ 可能的flag1 flag2 flag3)
                # 6字段: id resid type x y z
                # 9字段: id resid type x y z flag1 flag2 flag3
                sd = {"id": 0, "resid": 1, "type": 2, "coord": 3}
            elif n in (7, 10):
                # full format: id resid type charge x y z (+ 可能的flag1 flag2 flag3)
                # 7字段: id resid type charge x y z
                # 10字段: id resid type charge x y z flag1 flag2 flag3
                sd = {"id": 0, "resid": 1, "type": 2, "charge": 3, "coord": 4}
            else:
                # fallback for other formats
                sd = {"id": 0, "resid": 1, "type": 2, "coord": 3}
        else:
            sd = self.style_dict

        has_charge = "charge" in sd
        has_resid = "resid" in sd

        # atom ids aren't necessarily sequential
        atom_ids = np.zeros(n_atoms, dtype=np.int32)
        types = np.zeros(n_atoms, dtype=object)
        coords = np.zeros((n_atoms, 3), dtype=np.float32)
        if has_resid:
            resids = np.zeros(n_atoms, dtype=np.int32)
        else:
            resids = np.ones(n_atoms, dtype=np.int32)
        if has_charge:
            charges = np.zeros(n_atoms, dtype=np.float32)

        for i, line in enumerate(datalines):
            line = line.split()

            # these numpy array are already typed correctly,
            # so just pass the raw strings
            # and let numpy handle the conversion
            atom_ids[i] = line[sd["id"]]
            if has_resid:
                resids[i] = line[sd["resid"]]
            types[i] = line[sd["type"]]
            if has_charge:
                charges[i] = line[sd["charge"]]
            if "coord" in sd:
                # 自动检测格式使用coord字段
                coords[i] = np.array(
                    list(map(np.float32, line[sd["coord"] : sd["coord"] + 3]))
                )
            else:
                # 明确指定格式使用x, y, z字段
                coords[i] = np.array(
                    [
                        np.float32(line[sd["x"]]),
                        np.float32(line[sd["y"]]),
                        np.float32(line[sd["z"]]),
                    ]
                )
        # at this point, we've read the atoms section,
        # but it's still (potentially) unordered
        # TODO: Maybe we can optimise by checking if we need to sort
        # ie `if np.any(np.diff(atom_ids) > 1)`  but we want to search
        # in a generatorish way, np.any() would check everything at once
        order = np.argsort(atom_ids)
        atom_ids = atom_ids[order]
        types = types[order]
        coords = coords[order]
        if has_resid:
            resids = resids[order]
        if has_charge:
            charges = charges[order]

        atm = "Atom_Base"
        mle = "Molecule_Base"
        dbase = Database(n_atoms)
        Coordinate(coords, sid=atm, database=dbase)
        species = Species(types, sid=atm, database=dbase)
        if has_charge:
            Charge(charges, sid=atm, database=dbase)
        masses = np.zeros(n_atoms, dtype=np.float32)
        for i, at in enumerate(types):
            masses[i] = massdict[at]
        Mass(masses, sid=atm, database=dbase)
        Element.start_from_masses(masses, sid=atm, database=dbase)
        ids = ID(atom_ids, sid=(atm,), database=dbase)
        residx, resids = squash_by(resids)[:2]
        ids._update_source(resids, (mle,))
        dbase.ix._update_source(np.arange(len(resids)), (mle,))

        n_residues = len(resids)
        compomtrx = self.residx2mtrx(residx)
        compo = Composition(
            compomtrx, sid=(mle, atm), database=dbase, N=n_residues, M=n_atoms
        )

        return ids, species, compo, dbase, order, masses

    @staticmethod
    def residx2mtrx(residx):
        M = len(residx)
        data = np.array([1 for _ in range(M)])
        row = residx
        col = np.array([_ for _ in range(M)])
        mtrx = np.array([row, col, data])
        return mtrx

    def _parse_masses(self, datalines):
        """Lammps defines mass on a per atom type basis.

        This reads mass for each type and stores in dict
        """
        masses = {}
        for line in datalines:
            line = line.split()
            masses[line[0]] = float(line[1])

        return masses

    def _parse_box(self, header):
        x1, x2 = np.float32(header["xlo xhi"].split())
        x = x2 - x1
        y1, y2 = np.float32(header["ylo yhi"].split())
        y = y2 - y1
        z1, z2 = np.float32(header["zlo zhi"].split())
        z = z2 - z1

        if "xy xz yz" in header:
            # Triclinic
            unitcell = np.zeros((3, 3), dtype=np.float32)

            xy, xz, yz = np.float32(header["xy xz yz"].split())

            unitcell[0][0] = x
            unitcell[1][0] = xy
            unitcell[1][1] = y
            unitcell[2][0] = xz
            unitcell[2][1] = yz
            unitcell[2][2] = z

            # unitcell = triclinic_box(*unitcell)
        else:
            # Orthogonal
            unitcell = np.zeros(12, dtype=np.float32)
            unitcell[:3] = x, y, z
            unitcell[3:6] = 90.0, 90.0, 90.0
            unitcell[6:] = x1, x2, y1, y2, z1, z2

        return unitcell


class LAMMPSReader(ReaderBase):
    pass
