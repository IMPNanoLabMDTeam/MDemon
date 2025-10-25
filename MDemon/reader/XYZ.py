"""
Extended XYZ Format Reader for MDemon

This reader supports the Extended XYZ format as used by GPUMD and other molecular dynamics codes.
Format specification: https://gpumd.org/gpumd/input_files/model_xyz.html

Extended XYZ Format:
--------------------
Line 1: Number of atoms (integer)
Line 2: Keyword=value pairs (Lattice, Properties, pbc, etc.)
Line 3+: Atom data according to Properties specification

Example:
--------
80
Lattice="24.21 0.0 0.0 0.0 6.17 0.0 0.0 0.0 5.88" Properties=species:S:1:pos:R:3
Ga 1.918494 1.542860 1.380208
O 6.004629 0.000000 0.047220
...
"""

import re

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
    Velocity,
)
from ..core.universe import Box
from .base import ReaderBase, squash_by


class XYZReader(ReaderBase):
    """
    Reader for Extended XYZ format files.

    Supports the following keywords in line 2:
    - Lattice="ax ay az bx by bz cx cy cz" (Mandatory)
    - Properties=prop1:type1:ncols1:prop2:type2:ncols2:... (Mandatory)
    - pbc="T/F T/F T/F" (Optional, default: "T T T")

    Supported properties:
    - species:S:1 - Atom type/species (Mandatory)
    - pos:R:3 - Position vector (Mandatory)
    - mass:R:1 - Mass (Optional)
    - charge:R:1 - Charge (Optional)
    - vel:R:3 - Velocity vector (Optional)
    - group:I:N - Grouping methods (Optional, N can be any positive integer)
    """

    format = ["XYZ", "EXYZ"]

    def __init__(self, filename):
        super().__init__(filename)
        self.n_atoms = 0
        self.lattice = None
        self.pbc = None
        self.property_spec = {}

    def parse(self, **kwargs):
        """
        Parse Extended XYZ file and create Database.

        Returns
        -------
        database : Database
            The populated database object
        """
        with open(self.filename, "r") as f:
            lines = f.readlines()

        # Parse line 1: number of atoms
        self.n_atoms = int(lines[0].strip())

        # Parse line 2: metadata (Lattice, Properties, pbc, etc.)
        metadata = self._parse_metadata(lines[1].strip())

        # Parse atomic data (line 3 onwards)
        atom_data = self._parse_atom_data(lines[2 : 2 + self.n_atoms])

        # Create Database and populate attributes
        dbase = self._create_database(atom_data)

        return dbase

    def _parse_metadata(self, line):
        """
        Parse the metadata line (line 2) containing Lattice, Properties, etc.

        Parameters
        ----------
        line : str
            The metadata line

        Returns
        -------
        metadata : dict
            Dictionary containing parsed metadata
        """
        metadata = {}

        # Regular expression to match keyword="value" or keyword=value
        # Handles both quoted and unquoted values
        pattern = r'(\w+)\s*=\s*"([^"]*)"|(\w+)\s*=\s*(\S+)'
        matches = re.finditer(pattern, line, re.IGNORECASE)

        for match in matches:
            if match.group(1):  # Quoted value
                key = match.group(1).lower()
                value = match.group(2)
            else:  # Unquoted value
                key = match.group(3).lower()
                value = match.group(4)
            metadata[key] = value

        # Parse Lattice (Mandatory)
        if "lattice" not in metadata:
            raise ValueError("Lattice keyword is mandatory in Extended XYZ format")
        self.lattice = self._parse_lattice(metadata["lattice"])

        # Parse Properties (Mandatory)
        if "properties" not in metadata:
            raise ValueError("Properties keyword is mandatory in Extended XYZ format")
        self.property_spec = self._parse_properties(metadata["properties"])

        # Parse pbc (Optional, default: T T T)
        if "pbc" in metadata:
            pbc_str = metadata["pbc"].upper().split()
            self.pbc = [p == "T" for p in pbc_str]
        else:
            self.pbc = [True, True, True]

        return metadata

    def _parse_lattice(self, lattice_str):
        """
        Parse lattice string into 3x3 matrix.

        Parameters
        ----------
        lattice_str : str
            String containing 9 numbers: "ax ay az bx by bz cx cy cz"

        Returns
        -------
        lattice : np.ndarray
            3x3 lattice matrix where rows are [a, b, c] vectors
        """
        values = list(map(float, lattice_str.split()))
        if len(values) != 9:
            raise ValueError(f"Lattice must have 9 values, got {len(values)}")

        # Reshape to 3x3 matrix (row-major: [a, b, c])
        lattice = np.array(values, dtype=np.float32).reshape(3, 3)
        return lattice

    def _parse_properties(self, properties_str):
        """
        Parse Properties string to determine column structure.

        Format: prop1:type1:ncols1:prop2:type2:ncols2:...

        Parameters
        ----------
        properties_str : str
            Properties specification string

        Returns
        -------
        property_spec : dict
            Dictionary mapping property name to (data_type, n_columns, start_col)
            Example: {'species': ('S', 1, 0), 'pos': ('R', 3, 1)}
        """
        property_spec = {}
        parts = properties_str.split(":")

        if len(parts) % 3 != 0:
            raise ValueError(f"Properties format error: {properties_str}")

        col_index = 0
        for i in range(0, len(parts), 3):
            prop_name = parts[i].lower()
            data_type = parts[i + 1].upper()  # S=string, R=real, I=integer
            n_columns = int(parts[i + 2])

            property_spec[prop_name] = {
                "type": data_type,
                "ncols": n_columns,
                "start": col_index,
            }
            col_index += n_columns

        # Validate mandatory properties
        if "species" not in property_spec:
            raise ValueError("'species' property is mandatory")
        if "pos" not in property_spec:
            raise ValueError("'pos' property is mandatory")

        return property_spec

    def _parse_atom_data(self, lines):
        """
        Parse atomic data lines according to property specification.

        Parameters
        ----------
        lines : list of str
            Atomic data lines

        Returns
        -------
        atom_data : dict
            Dictionary containing arrays for each property
        """
        if len(lines) < self.n_atoms:
            raise ValueError(f"Expected {self.n_atoms} atom lines, found {len(lines)}")

        # Initialize data arrays
        atom_data = {
            "species": [],
            "pos": np.zeros((self.n_atoms, 3), dtype=np.float32),
        }

        # Optional properties
        if "mass" in self.property_spec:
            atom_data["mass"] = np.zeros(self.n_atoms, dtype=np.float32)
        if "charge" in self.property_spec:
            atom_data["charge"] = np.zeros(self.n_atoms, dtype=np.float32)
        if "vel" in self.property_spec:
            atom_data["vel"] = np.zeros((self.n_atoms, 3), dtype=np.float32)
        if "group" in self.property_spec:
            n_groups = self.property_spec["group"]["ncols"]
            atom_data["group"] = np.zeros((self.n_atoms, n_groups), dtype=np.int32)

        # Parse each line
        for i, line in enumerate(lines[: self.n_atoms]):
            tokens = line.strip().split()

            # Species
            spec = self.property_spec["species"]
            atom_data["species"].append(tokens[spec["start"]])

            # Position
            spec = self.property_spec["pos"]
            atom_data["pos"][i] = [float(tokens[spec["start"] + j]) for j in range(3)]

            # Optional: Mass
            if "mass" in self.property_spec:
                spec = self.property_spec["mass"]
                atom_data["mass"][i] = float(tokens[spec["start"]])

            # Optional: Charge
            if "charge" in self.property_spec:
                spec = self.property_spec["charge"]
                atom_data["charge"][i] = float(tokens[spec["start"]])

            # Optional: Velocity
            if "vel" in self.property_spec:
                spec = self.property_spec["vel"]
                atom_data["vel"][i] = [
                    float(tokens[spec["start"] + j]) for j in range(3)
                ]

            # Optional: Group
            if "group" in self.property_spec:
                spec = self.property_spec["group"]
                for j in range(spec["ncols"]):
                    atom_data["group"][i, j] = int(tokens[spec["start"] + j])

        return atom_data

    def _create_database(self, atom_data):
        """
        Create Database and populate with parsed data.

        Parameters
        ----------
        atom_data : dict
            Dictionary containing parsed atomic data

        Returns
        -------
        database : Database
            Populated database object
        """
        n_atoms = self.n_atoms
        dbase = Database(n_atoms)

        atm = "Atom_Base"
        mle = "Molecule_Base"

        # 1. Coordinates
        Coordinate(atom_data["pos"], sid=atm, database=dbase)

        # 2. Species - convert element symbols to integer types
        species_symbols = atom_data["species"]
        unique_species = []
        species_map = {}
        species_types = np.zeros(n_atoms, dtype=np.int32)

        for i, symbol in enumerate(species_symbols):
            if symbol not in species_map:
                species_map[symbol] = len(unique_species) + 1
                unique_species.append(symbol)
            species_types[i] = species_map[symbol]

        # Store species as integers (1-indexed)
        species = Species(species_types, sid=atm, database=dbase)

        # 3. Mass - use provided values or default from periodic table
        if "mass" in atom_data:
            masses = atom_data["mass"]
        else:
            # Get default masses from element symbols
            import periodictable

            masses = np.zeros(n_atoms, dtype=np.float32)
            for i, symbol in enumerate(species_symbols):
                try:
                    element = getattr(periodictable, symbol)
                    masses[i] = element.mass
                except (AttributeError, TypeError):
                    # If element not found, use species type as mass (placeholder)
                    masses[i] = species_types[i] * 12.0

        Mass(masses, sid=atm, database=dbase)

        # 4. Element - derive from masses
        Element.start_from_masses(masses, sid=atm, database=dbase)

        # 5. Charge (optional)
        if "charge" in atom_data:
            Charge(atom_data["charge"], sid=atm, database=dbase)

        # 6. Velocity (optional)
        if "vel" in atom_data:
            Velocity(atom_data["vel"], sid=atm, database=dbase)
        else:
            # Create zero velocities
            velocities = np.zeros((n_atoms, 3), dtype=np.float32)
            Velocity(velocities, sid=atm, database=dbase)

        # 7. IDs (sequential)
        atom_ids = np.arange(1, n_atoms + 1, dtype=np.int32)
        ids = ID(atom_ids, sid=(atm,), database=dbase)

        # 8. Create molecule structure (all atoms in one molecule by default)
        # This can be modified if group information suggests molecular structure
        resids = np.ones(n_atoms, dtype=np.int32)
        residx, resids_unique = squash_by(resids)[:2]
        ids._update_source(resids_unique, (mle,))
        dbase.ix._update_source(np.arange(len(resids_unique)), (mle,))

        # 9. Composition matrix (molecule-atom relationship)
        n_residues = len(resids_unique)
        row = residx
        col = np.arange(n_atoms)
        data = np.ones(n_atoms)
        compomtrx = np.array([row, col, data])

        compo = Composition(
            compomtrx, sid=(mle, atm), database=dbase, N=n_residues, M=n_atoms
        )

        # 10. Box information
        box = self._create_box_from_lattice()
        Box(box, database=dbase)

        return dbase

    def _create_box_from_lattice(self):
        """
        Convert lattice matrix to box format.

        For orthogonal boxes: [lx, ly, lz, 90, 90, 90, xlo, xhi, ylo, yhi, zlo, zhi]
        For triclinic boxes: 3x3 matrix

        Returns
        -------
        box : np.ndarray
            Box specification
        """
        # Check if lattice is orthogonal
        is_orthogonal = (
            np.abs(self.lattice[0, 1]) < 1e-6
            and np.abs(self.lattice[0, 2]) < 1e-6
            and np.abs(self.lattice[1, 0]) < 1e-6
            and np.abs(self.lattice[1, 2]) < 1e-6
            and np.abs(self.lattice[2, 0]) < 1e-6
            and np.abs(self.lattice[2, 1]) < 1e-6
        )

        if is_orthogonal:
            # Orthogonal box
            lx = self.lattice[0, 0]
            ly = self.lattice[1, 1]
            lz = self.lattice[2, 2]

            box = np.zeros(12, dtype=np.float32)
            box[:3] = [lx, ly, lz]
            box[3:6] = [90.0, 90.0, 90.0]  # angles
            box[6:] = [0.0, lx, 0.0, ly, 0.0, lz]  # lo, hi for each dimension
        else:
            # Triclinic box - use lattice matrix directly
            box = self.lattice.copy()

        return box


class XYZWriter:
    """
    Writer for Extended XYZ format files.

    This will write MDemon Universe data to Extended XYZ format.
    """

    def __init__(self, filename):
        self.filename = filename

    def write(self, universe, properties=None, include_velocities=False):
        """
        Write Universe to Extended XYZ file.

        Parameters
        ----------
        universe : Universe
            MDemon Universe object
        properties : list, optional
            List of properties to include. Default: ['species', 'pos']
        include_velocities : bool, optional
            Whether to include velocities. Default: False
        """
        if properties is None:
            properties = ["species", "pos"]
            if include_velocities:
                properties.append("vel")

        with open(self.filename, "w") as f:
            # Line 1: Number of atoms
            n_atoms = len(universe.atoms)
            f.write(f"{n_atoms}\n")

            # Line 2: Lattice and Properties
            box = universe.box
            if len(box.shape) == 1 and len(box) == 12:
                # Orthogonal box
                lattice_str = (
                    f"{box[0]:.10f} 0.0 0.0 "
                    f"0.0 {box[1]:.10f} 0.0 "
                    f"0.0 0.0 {box[2]:.10f}"
                )
            else:
                # Triclinic box
                lattice_str = " ".join([f"{v:.10f}" for v in box.flatten()])

            # Build properties string
            prop_specs = []
            if "species" in properties:
                prop_specs.append("species:S:1")
            if "pos" in properties:
                prop_specs.append("pos:R:3")
            if "vel" in properties:
                prop_specs.append("vel:R:3")

            properties_str = ":".join(prop_specs)

            f.write(f'Lattice="{lattice_str}" Properties={properties_str}\n')

            # Lines 3+: Atomic data
            for atom in universe.atoms:
                line_parts = []

                # Species (use element symbol or species number)
                if "species" in properties:
                    # Try to get element symbol, fallback to species number
                    try:
                        import periodictable

                        element = periodictable.elements[atom.element]
                        line_parts.append(str(element.symbol))
                    except:
                        line_parts.append(f"Type{atom.species}")

                # Position
                if "pos" in properties:
                    pos = atom.coordinate
                    line_parts.extend(
                        [f"{pos[0]:.10f}", f"{pos[1]:.10f}", f"{pos[2]:.10f}"]
                    )

                # Velocity
                if "vel" in properties:
                    vel = atom.velocity
                    line_parts.extend(
                        [f"{vel[0]:.10f}", f"{vel[1]:.10f}", f"{vel[2]:.10f}"]
                    )

                f.write(" ".join(line_parts) + "\n")
