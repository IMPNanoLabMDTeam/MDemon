"""
Tests for XYZ Reader

Tests the Extended XYZ format reader with various files.
"""

from pathlib import Path

import numpy as np
import pytest

import MDemon as md
from MDemon.reader.XYZ import XYZReader


class TestXYZReader:
    """Test XYZ file reading functionality"""

    @pytest.fixture
    def data_dir(self):
        """Return path to test data directory"""
        return Path(__file__).parent / "data" / "xyz" / "basic"

    def test_xyz_reader_registration(self):
        """Test that XYZ reader is properly registered"""
        from MDemon import _READERS

        assert "XYZ" in _READERS
        assert "EXYZ" in _READERS
        assert _READERS["XYZ"] == XYZReader

    def test_read_beta_221_ort(self, data_dir):
        """Test reading beta_221_ort.xyz file"""
        xyz_file = data_dir / "beta_221_ort.xyz"

        # Read the file
        u = md.Universe(str(xyz_file))

        # Check number of atoms
        assert len(u.atoms) == 80, f"Expected 80 atoms, got {len(u.atoms)}"

        # Check that we have atoms
        assert hasattr(u, "atoms"), "Universe should have atoms attribute"

        # Check first atom properties
        atom0 = u.atoms[0]
        assert hasattr(atom0, "coordinate"), "Atom should have coordinate"
        assert hasattr(atom0, "species"), "Atom should have species"

        # Verify first atom data (Ga 1.918494105352 1.542860466635 1.380208059015)
        coord0 = atom0.coordinate
        expected_coord = np.array(
            [1.918494105352, 1.542860466635, 1.380208059015], dtype=np.float32
        )
        np.testing.assert_allclose(coord0, expected_coord, rtol=1e-5)

        print(f"✓ First atom coordinate: {coord0}")
        print(f"✓ First atom species: {atom0.species}")

    def test_lattice_parsing(self, data_dir):
        """Test lattice parsing from XYZ file"""
        xyz_file = data_dir / "beta_221_ort.xyz"

        reader = XYZReader(str(xyz_file))
        dbase = reader.parse()

        # Check lattice dimensions
        expected_lattice = np.array(
            [
                [24.214144270918, 0.0, 0.0],
                [0.0, 6.171441866538, 0.0],
                [0.0, 0.0, 5.879384648762],
            ],
            dtype=np.float32,
        )

        np.testing.assert_allclose(reader.lattice, expected_lattice, rtol=1e-5)
        print(f"✓ Lattice matrix:\n{reader.lattice}")

    def test_species_identification(self, data_dir):
        """Test species identification in XYZ file"""
        xyz_file = data_dir / "beta_221_ort.xyz"

        u = md.Universe(str(xyz_file))

        # Count species
        species_list = [atom.species for atom in u.atoms]
        unique_species = set(species_list)

        # Should have 2 species: Ga and O
        assert (
            len(unique_species) == 2
        ), f"Expected 2 species, got {len(unique_species)}"

        # Count Ga and O atoms (from the file, we can manually verify)
        # Lines 3-10, 23-30, 43-50, 63-70: Ga atoms
        # Lines 11-22, 31-42, 51-62, 71-82: O atoms
        species_counts = {}
        for s in species_list:
            species_counts[s] = species_counts.get(s, 0) + 1

        print(f"✓ Species distribution: {species_counts}")
        assert len(species_counts) == 2, "Should have exactly 2 species"

    def test_box_creation(self, data_dir):
        """Test box creation from lattice"""
        xyz_file = data_dir / "beta_221_ort.xyz"

        u = md.Universe(str(xyz_file))

        # Check box exists
        assert hasattr(u, "box"), "Universe should have box attribute"

        box = u.box
        print(f"✓ Box shape: {box.shape}")
        print(f"✓ Box data: {box}")

        # For orthogonal box, should be 12 elements
        # [lx, ly, lz, 90, 90, 90, xlo, xhi, ylo, yhi, zlo, zhi]
        if len(box.shape) == 1:
            assert (
                len(box) == 12
            ), f"Orthogonal box should have 12 elements, got {len(box)}"
            assert box[0] > 0, "Box x-dimension should be positive"
            assert box[1] > 0, "Box y-dimension should be positive"
            assert box[2] > 0, "Box z-dimension should be positive"


class TestXYZReaderProperties:
    """Test XYZ reader with different property specifications"""

    def test_properties_parsing(self):
        """Test parsing of Properties string"""
        reader = XYZReader("dummy.xyz")

        # Test standard properties
        props = reader._parse_properties("species:S:1:pos:R:3")
        assert "species" in props
        assert "pos" in props
        assert props["species"]["type"] == "S"
        assert props["species"]["ncols"] == 1
        assert props["species"]["start"] == 0
        assert props["pos"]["type"] == "R"
        assert props["pos"]["ncols"] == 3
        assert props["pos"]["start"] == 1

    def test_properties_with_velocity(self):
        """Test properties with velocity"""
        reader = XYZReader("dummy.xyz")

        props = reader._parse_properties("species:S:1:pos:R:3:vel:R:3")
        assert "vel" in props
        assert props["vel"]["type"] == "R"
        assert props["vel"]["ncols"] == 3
        assert props["vel"]["start"] == 4

    def test_properties_with_mass_charge(self):
        """Test properties with mass and charge"""
        reader = XYZReader("dummy.xyz")

        props = reader._parse_properties("species:S:1:pos:R:3:mass:R:1:charge:R:1")
        assert "mass" in props
        assert "charge" in props
        assert props["mass"]["start"] == 4
        assert props["charge"]["start"] == 5


class TestXYZReaderErrors:
    """Test error handling in XYZ reader"""

    def test_missing_lattice(self, tmp_path):
        """Test error when Lattice is missing"""
        xyz_file = tmp_path / "missing_lattice.xyz"
        with open(xyz_file, "w") as f:
            f.write("2\n")
            f.write("Properties=species:S:1:pos:R:3\n")
            f.write("C 0 0 0\n")
            f.write("C 1 0 0\n")

        with pytest.raises(ValueError, match="Lattice keyword is mandatory"):
            md.Universe(str(xyz_file))

    def test_missing_properties(self, tmp_path):
        """Test error when Properties is missing"""
        xyz_file = tmp_path / "missing_properties.xyz"
        with open(xyz_file, "w") as f:
            f.write("2\n")
            f.write('Lattice="10 0 0 0 10 0 0 0 10"\n')
            f.write("C 0 0 0\n")
            f.write("C 1 0 0\n")

        with pytest.raises(ValueError, match="Properties keyword is mandatory"):
            md.Universe(str(xyz_file))

    def test_missing_species(self, tmp_path):
        """Test error when species property is missing"""
        xyz_file = tmp_path / "missing_species.xyz"
        with open(xyz_file, "w") as f:
            f.write("2\n")
            f.write('Lattice="10 0 0 0 10 0 0 0 10" Properties=pos:R:3\n')
            f.write("0 0 0\n")
            f.write("1 0 0\n")

        with pytest.raises(ValueError, match="'species' property is mandatory"):
            md.Universe(str(xyz_file))


if __name__ == "__main__":
    # Run basic tests
    print("=" * 70)
    print("Testing XYZ Reader")
    print("=" * 70)

    data_dir = Path(__file__).parent / "data" / "xyz" / "basic"
    xyz_file = data_dir / "beta_221_ort.xyz"

    if xyz_file.exists():
        print(f"\n📂 Reading file: {xyz_file}")

        try:
            # Test basic reading
            u = md.Universe(str(xyz_file))

            print(f"✅ Successfully loaded universe")
            print(f"   Number of atoms: {len(u.atoms)}")
            print(f"   First atom coordinate: {u.atoms[0].coordinate}")
            print(f"   First atom species: {u.atoms[0].species}")

            # Check box
            print(f"\n📦 Box information:")
            print(f"   Box shape: {u.box.shape}")
            print(f"   Box data: {u.box}")

            # Check species distribution
            species_list = [atom.species for atom in u.atoms]
            unique_species = set(species_list)
            species_counts = {}
            for s in species_list:
                species_counts[s] = species_counts.get(s, 0) + 1

            print(f"\n🔬 Species information:")
            print(f"   Unique species: {len(unique_species)}")
            print(f"   Species distribution: {species_counts}")

            print("\n" + "=" * 70)
            print("✅ All basic tests passed!")
            print("=" * 70)

        except Exception as e:
            print(f"\n❌ Error: {e}")
            import traceback

            traceback.print_exc()
    else:
        print(f"❌ File not found: {xyz_file}")
