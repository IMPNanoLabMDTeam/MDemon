"""
Test suite for MDemon constants module.

This module tests the physical constants, unit conversions, and utility functions.
"""

import pytest
import math
from MDemon.constants import PHYSICS, CHEMISTRY, CONVERSION, MATERIAL
from MDemon.constants import k_B, N_A, c_light, ANGSTROM_TO_METER, EV_TO_JOULE
from MDemon.constants import utils


class TestPhysicsConstants:
    """Test fundamental physical constants"""
    
    def test_boltzmann_constant(self):
        """Test Boltzmann constant value"""
        assert abs(k_B - 1.380649e-23) < 1e-28
        assert abs(PHYSICS.k_B - 1.380649e-23) < 1e-28
    
    def test_avogadro_constant(self):
        """Test Avogadro's number"""
        assert abs(N_A - 6.02214076e23) < 1e16
        assert abs(PHYSICS.N_A - 6.02214076e23) < 1e16
    
    def test_speed_of_light(self):
        """Test speed of light"""
        assert c_light == 299792458.0
        assert PHYSICS.c == 299792458.0
    
    def test_elementary_charge(self):
        """Test elementary charge"""
        assert abs(PHYSICS.e - 1.602176634e-19) < 1e-28


class TestChemistryConstants:
    """Test chemistry-related constants"""
    
    def test_gas_constant(self):
        """Test gas constant consistency"""
        # R = k_B * N_A
        calculated_R = PHYSICS.k_B * PHYSICS.N_A
        assert abs(CHEMISTRY.R - calculated_R) < 1e-10
    
    def test_faraday_constant(self):
        """Test Faraday constant"""
        # F = e * N_A
        calculated_F = PHYSICS.e * PHYSICS.N_A
        assert abs(CHEMISTRY.F - calculated_F) < 1e-5
    
    def test_bond_lengths(self):
        """Test bond length lookups"""
        cc_bond = CHEMISTRY.get_bond_length('C-C')
        assert cc_bond is not None
        assert 1.5e-10 < cc_bond < 1.6e-10  # Should be around 1.54 Å
        
        # Test non-existent bond
        unknown_bond = CHEMISTRY.get_bond_length('X-Y')
        assert unknown_bond is None
    
    def test_vdw_radii(self):
        """Test van der Waals radii"""
        h_radius = CHEMISTRY.get_vdw_radius('H')
        assert h_radius is not None
        assert 1.0e-10 < h_radius < 1.3e-10  # Should be around 1.2 Å
        
        # Test non-existent element
        unknown_radius = CHEMISTRY.get_vdw_radius('Xx')
        assert unknown_radius is None


class TestConversions:
    """Test unit conversion factors"""
    
    def test_length_conversions(self):
        """Test length conversion factors"""
        assert CONVERSION.ANGSTROM_TO_METER == 1e-10
        assert CONVERSION.NANOMETER_TO_METER == 1e-9
        assert CONVERSION.METER_TO_ANGSTROM == 1e10
    
    def test_energy_conversions(self):
        """Test energy conversion factors"""
        # Test eV to Joule
        expected_ev_to_j = 1.602176634e-19
        assert abs(CONVERSION.EV_TO_JOULE - expected_ev_to_j) < 1e-28
        
        # Test Hartree to eV
        expected_ha_to_ev = 27.211386245988
        assert abs(CONVERSION.HARTREE_TO_EV - expected_ha_to_ev) < 1e-10
    
    def test_temperature_conversions(self):
        """Test temperature conversion functions"""
        # 0°C should be 273.15 K
        temp_k = CONVERSION.temperature_celsius_to_kelvin(0.0)
        assert abs(temp_k - 273.15) < 1e-10
        
        # 32°F should be 273.15 K
        temp_k_f = CONVERSION.temperature_fahrenheit_to_kelvin(32.0)
        assert abs(temp_k_f - 273.15) < 1e-10


class TestMaterialProperties:
    """Test material properties"""
    
    def test_ga2o3_properties(self):
        """Test Ga₂O₃ properties"""
        # Test β-Ga₂O₃ density
        density = MATERIAL.get_density('Ga2O3_beta')
        assert density is not None
        assert 5000 < density < 6500  # Should be around 5950 kg/m³
        
        # Test lattice parameters
        params = MATERIAL.get_lattice_parameters('Ga2O3_beta')
        assert params is not None
        assert 'a' in params and 'b' in params and 'c' in params and 'beta' in params
        assert 12 < params['a'] < 13  # Around 12.214 Å
        assert 3 < params['b'] < 4    # Around 3.037 Å
    
    def test_common_materials(self):
        """Test common material properties"""
        # Silicon properties
        si_density = MATERIAL.get_density('Si')
        assert si_density is not None
        assert 2000 < si_density < 2500
        
        # Copper properties
        cu_thermal_cond = MATERIAL.get_thermal_conductivity('Cu')
        assert cu_thermal_cond is not None
        assert 300 < cu_thermal_cond < 450
    
    def test_material_listing(self):
        """Test material listing function"""
        materials = MATERIAL.list_available_materials()
        assert isinstance(materials, list)
        assert 'Si' in materials
        assert 'Cu' in materials
        assert 'Ga2O3_beta' in materials


class TestUtilityFunctions:
    """Test utility functions"""
    
    def test_energy_conversion_utils(self):
        """Test energy conversion utility functions"""
        # Convert 1 eV to Joules
        energy_j = utils.convert_energy(1.0, 'eV', 'J')
        expected = CONVERSION.EV_TO_JOULE
        assert abs(energy_j - expected) < 1e-25
        
        # Convert 1 Hartree to eV
        energy_ev = utils.convert_energy(1.0, 'Hartree', 'eV')
        expected = CONVERSION.HARTREE_TO_EV
        assert abs(energy_ev - expected) < 1e-10
        
        # Test round-trip conversion
        original = 2.5  # eV
        converted = utils.convert_energy(original, 'eV', 'kcal/mol')
        back = utils.convert_energy(converted, 'kcal/mol', 'eV')
        assert abs(back - original) < 1e-10
    
    def test_length_conversion_utils(self):
        """Test length conversion utility functions"""
        # Convert 5 Å to nanometers
        length_nm = utils.convert_length(5.0, 'Å', 'nm')
        assert abs(length_nm - 0.5) < 1e-10
        
        # Convert 1 bohr to Ångström
        length_ang = utils.convert_length(1.0, 'bohr', 'Å')
        expected = CONVERSION.BOHR_TO_METER / CONVERSION.ANGSTROM_TO_METER
        assert abs(length_ang - expected) < 1e-10
    
    def test_temperature_conversion_utils(self):
        """Test temperature conversion utility functions"""
        # Room temperature conversions
        temp_c = 25.0
        temp_k = utils.convert_temperature(temp_c, 'C', 'K')
        temp_f = utils.convert_temperature(temp_c, 'C', 'F')
        
        assert abs(temp_k - 298.15) < 1e-10
        assert abs(temp_f - 77.0) < 1e-10
        
        # Test round-trip
        back_c = utils.convert_temperature(temp_k, 'K', 'C')
        assert abs(back_c - temp_c) < 1e-10
    
    def test_boltzmann_energy_calculation(self):
        """Test thermal energy calculation"""
        temp = 298.15  # K
        
        # Calculate k_B * T in different units
        energy_j = utils.boltzmann_energy(temp, 'J')
        energy_ev = utils.boltzmann_energy(temp, 'eV')
        
        expected_j = PHYSICS.k_B * temp
        expected_ev = expected_j / CONVERSION.EV_TO_JOULE
        
        assert abs(energy_j - expected_j) < 1e-25
        assert abs(energy_ev - expected_ev) < 1e-15
    
    def test_thermal_velocity(self):
        """Test thermal velocity calculation"""
        # Hydrogen atom at room temperature
        mass_h = 1.008 * CHEMISTRY.u  # kg
        temp = 298.15  # K
        
        v_thermal = utils.thermal_velocity(mass_h, temp)
        
        # Expected: v = sqrt(3 * k_B * T / m)
        expected = math.sqrt(3 * PHYSICS.k_B * temp / mass_h)
        assert abs(v_thermal - expected) < 1e-5
        
        # Should be reasonable velocity for hydrogen (~2500 m/s)
        assert 2000 < v_thermal < 3000
    
    def test_planck_frequency(self):
        """Test frequency calculation from energy"""
        # Visible light photon (2 eV)
        energy_j = 2.0 * CONVERSION.EV_TO_JOULE
        frequency = utils.planck_frequency(energy_j, 'Hz')
        
        expected = energy_j / PHYSICS.h
        assert abs(frequency - expected) < 1e5  # Allow some numerical error
        
        # Should be in visible light range (~10^15 Hz)
        assert 1e14 < frequency < 1e16
    
    def test_material_property_lookup(self):
        """Test material property lookup utility"""
        # Test valid material and property
        density = utils.get_material_property('Si', 'density')
        assert density is not None
        assert 2000 < density < 2500
        
        # Test invalid property
        with pytest.raises(ValueError):
            utils.get_material_property('Si', 'invalid_property')
        
        # Test non-existent material
        result = utils.get_material_property('NonExistentMaterial', 'density')
        assert result is None
    
    def test_debye_length(self):
        """Test Debye screening length calculation"""
        temp = 298.15  # K
        ionic_strength = 0.1  # mol/L
        
        debye_len = utils.debye_length(temp, ionic_strength, 'm')
        
        # Should be reasonable for physiological conditions (~1 nm)
        assert 5e-10 < debye_len < 2e-9
        
        # Test unit conversion
        debye_len_nm = utils.debye_length(temp, ionic_strength, 'nm')
        expected_nm = debye_len / 1e-9
        assert abs(debye_len_nm - expected_nm) < 1e-12
    
    def test_available_units(self):
        """Test available units listing"""
        units = utils.list_available_units()
        assert isinstance(units, dict)
        assert 'energy' in units
        assert 'length' in units
        assert 'temperature' in units
        
        assert 'eV' in units['energy']
        assert 'Å' in units['length']
        assert 'K' in units['temperature']


class TestIntegration:
    """Integration tests for the constants system"""
    
    def test_consistency_checks(self):
        """Test consistency between related constants"""
        # Gas constant should equal k_B * N_A
        calculated_R = PHYSICS.k_B * PHYSICS.N_A
        assert abs(CHEMISTRY.R - calculated_R) < 1e-10
        
        # Faraday constant should equal e * N_A
        calculated_F = PHYSICS.e * PHYSICS.N_A
        assert abs(CHEMISTRY.F - calculated_F) < 1e-5
        
        # Reduced Planck constant
        calculated_hbar = PHYSICS.h / (2 * math.pi)
        assert abs(PHYSICS.hbar - calculated_hbar) < 1e-40
    
    def test_unit_consistency(self):
        """Test that conversions are consistent"""
        # Test energy conversion consistency
        energy_ha = 1.0  # Hartree
        energy_ev = utils.convert_energy(energy_ha, 'Hartree', 'eV')
        energy_j = utils.convert_energy(energy_ha, 'Hartree', 'J')
        
        # Convert eV to J and compare
        energy_j_via_ev = utils.convert_energy(energy_ev, 'eV', 'J')
        assert abs(energy_j - energy_j_via_ev) < 1e-25
        
        # Test length conversion consistency
        length_bohr = 1.0
        length_ang = utils.convert_length(length_bohr, 'bohr', 'Å')
        length_m = utils.convert_length(length_bohr, 'bohr', 'm')
        
        # Convert Å to m and compare
        length_m_via_ang = utils.convert_length(length_ang, 'Å', 'm')
        assert abs(length_m - length_m_via_ang) < 1e-20


if __name__ == "__main__":
    # Run tests if script is executed directly
    pytest.main([__file__, "-v"]) 