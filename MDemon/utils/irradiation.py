import math
import numpy as np
import json
from datetime import datetime
import sys

from ..constants import PHYSICS, IRRADIATION, CONVERSION, MATERIAL


class WaligorskiZhangCalculator:
    """
    Waligorski-Zhang radial dose distribution model calculator
    
    This calculator uses constants from the MDemon.constants module and
    requires material properties to be specified for calculations.
    """
    
    def __init__(self, molecular_weight, atoms_per_molecule, density_g_per_cm3,ion_name, ion_energy_MeV_per_amu, ion_Z, energy_loss_keV_per_um,
                 g_factor = 0.17):
        """
        Initialize the calculator with material properties.
        
        Parameters:
        -----------
        molecular_weight : float
            Molecular weight of the target material (for example 18 g/mol for H2O)
        atoms_per_molecule : int
            Number of atoms per molecule in the target material (for example 3 for H2O)
        density_g_per_cm3 : float
            Density of the target material (for example 1.0 for H2O)
        ion_name : str
            Name of the ion (for example 'Ga' for Ga2O3)
        ion_energy_MeV_per_amu : float
            Energy of the ion per atomic mass unit (for example 1 MeV/amu for Xe ion)
        ion_Z : int
            Atomic number of the ion (for example 54 for Xe)
        energy_loss_keV_per_um : float
            Energy loss per micrometer of the ion in the target material (for example 10 keV/um for Xe in Ga2O3)
        g_factor : float
            g-factor of the ion (for example 0.17 for Xe in PET)
        """
        # Import constants from MDemon.constants
        self.electron_mass_keV = IRRADIATION.electron_mass_keV
        self.water_constant_keV_per_mm = IRRADIATION.water_constant_keV_per_mm
        self.avogadro_constant = PHYSICS.N_A
        self.eV_to_K = IRRADIATION.eV_to_K_conversion
        self.amu_to_MeV = IRRADIATION.amu_to_MeV
        self.barkas_coefficient = IRRADIATION.barkas_coefficient
        self.ionization_energy_eV = IRRADIATION.ionization_energy_eV
        self.k_g_per_cm2_keV_per_alpha = IRRADIATION.wz_range_g_per_cm2_keV_per_alpha
        
        # Target material properties
        self.molecular_weight = molecular_weight
        self.atoms_per_molecule = atoms_per_molecule
        self.density_g_per_cm3 = density_g_per_cm3

        # SHI(swift heavy ion) information
        self.ion_name = ion_name
        self.ion_energy_MeV_per_amu = ion_energy_MeV_per_amu
        self.ion_Z = ion_Z
        self.energy_loss_keV_per_um = energy_loss_keV_per_um
        
    def calculate_beta_gamma(self, energy_MeV_per_amu):
        """Calculate relativistic β and γ factors"""
        γ = 1 + energy_MeV_per_amu / self.amu_to_MeV
        β = math.sqrt(1 - 1/γ**2)
        return β, γ
    
    def calculate_effective_charge(self, Z, β):
        """Calculate effective charge using Barkas correction"""
        return Z * (1 - math.exp(-self.barkas_coefficient * β * Z**(-2/3)))
    
    def calculate_range_parameter_alpha(self, β):
        """Calculate range parameter α"""
        return IRRADIATION.get_range_parameter_alpha(β)
    
    def calculate_radial_dose(self, radius_nm):
        """Calculate radial dose distribution"""
        β, γ = self.calculate_beta_gamma(self.ion_energy_MeV_per_amu)
        Z_eff = self.calculate_effective_charge(self.ion_Z, β)
        α = self.calculate_range_parameter_alpha(β)
        
        ionization_energy_keV = self.ionization_energy_eV / 1000.0
        θ = self.k_g_per_cm2_keV_per_alpha * (ionization_energy_keV ** 1.079)
        W = 2 * self.electron_mass_keV * β**2 * γ**2
        T = self.k_g_per_cm2_keV_per_alpha * (W ** α)
        
        radius_cm = radius_nm * 1e-7
        t_g_per_cm2 = radius_cm * self.density_g_per_cm3
        
        if t_g_per_cm2 <= 0 or (T + θ) <= 0 or α <= 0:
            return 0
        
        water_constant_keV_per_cm = self.water_constant_keV_per_mm * 1e1
        constant_factor = water_constant_keV_per_cm * Z_eff**2 / (2 * np.pi * α * β**2 * t_g_per_cm2)
        fraction = (t_g_per_cm2 + θ) / (T + θ)
        
        if fraction >= 1:
            return 0
        
        power_term = (1 - fraction) ** (1/α)
        dose_keV_cm3_per_g2 = constant_factor * power_term / (t_g_per_cm2 + θ)
        dose_keV_per_cm3 = dose_keV_cm3_per_g2 * (self.density_g_per_cm3)**2
        dose_keV_per_nm3 = dose_keV_per_cm3 * 1e-21

        return dose_keV_per_nm3
    
    def calculate_radial_energy_distribution(self, radius_array_nm):
        """Calculate complete radial energy distribution"""
        dose_array = np.array([
            self.calculate_radial_dose(r)
            for r in radius_array_nm
        ])
        
        atomic_density_nm3 = (self.density_g_per_cm3 * self.avogadro_constant * 1e-21) / self.molecular_weight * self.atoms_per_molecule
        
        energy_per_atom_eV = (dose_array * 1000.0) / atomic_density_nm3
        
        # Calculate cumulative energy
        cumulative_energy = np.zeros_like(radius_array_nm)
        for i in range(1, len(radius_array_nm)):
            dr = radius_array_nm[i] - radius_array_nm[i-1]
            ring_area = 2 * np.pi * radius_array_nm[i] * dr
            ring_energy = dose_array[i] * ring_area * 1.0
            cumulative_energy[i] = cumulative_energy[i-1] + ring_energy

        # normalize the dose by the energy loss per um
        normalized_dose = dose_array * self.g_factor * self.energy_loss_keV_per_um / cumulative_energy[-1]

        # calculate the temperature
        temperature = energy_per_atom_eV * self.eV_to_K

        return {
            'radius': radius_array_nm,
            'dose_density': dose_array,
            'normalized_dose_density': normalized_dose,
            'energy_per_atom': energy_per_atom_eV,
            'temperature': temperature,
            'cumulative_energy': cumulative_energy,
            'atomic_density': atomic_density_nm3,
            'calculated_energy_loss': cumulative_energy[-1] if len(cumulative_energy) > 0 else 0
        }

def plot_results(figure, radius, dose, normalized_dose, temperature, g_factor, ion_name, material_name, analysis_info):
    """Plot calculation results"""
    figure.clear()
    
    ax1 = figure.add_subplot(221)
    ax2 = figure.add_subplot(222)
    ax3 = figure.add_subplot(223)
    ax4 = figure.add_subplot(224)
    
    # Subplot 1: Dose distribution - log scale
    ax1.loglog(radius, dose, 'b-', label='Original Dose')
    ax1.loglog(radius, normalized_dose, 'r--', label=f'Normalized (g={g_factor})')
    ax1.set_xlabel('Radial Distance (nm)')
    ax1.set_ylabel('Dose (keV/nm^3)')
    ax1.set_title('Dose Distribution (Log)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Subplot 2: Dose distribution - linear scale
    ax2.plot(radius, dose, 'b-', label='Original Dose')
    ax2.plot(radius, normalized_dose, 'r--', label=f'Normalized (g={g_factor})')
    ax2.set_xlabel('Radial Distance (nm)')
    ax2.set_ylabel('Dose (keV/nm^3)')
    ax2.set_title('Dose Distribution (Linear)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Subplot 3: Temperature distribution - log scale
    ax3.loglog(radius, temperature, 'g-', label='Temperature')
    ax3.set_xlabel('Radial Distance (nm)')
    ax3.set_ylabel('Temperature (K)')
    ax3.set_title('Temperature Distribution (Log)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # Subplot 4: Energy per atom distribution
    energy_per_atom = analysis_info['energy_per_atom_eV']
    ax4.loglog(radius, energy_per_atom, 'm-', label=f'Energy per {material_name} Atom')
    ax4.set_xlabel('Radial Distance (nm)')
    ax4.set_ylabel('Energy per Atom (eV)')
    ax4.set_title('Single Atom Energy Distribution')
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    main_title = f'{ion_name} Radiation in {material_name}'
    figure.suptitle(main_title, fontsize=14, fontweight='bold')
    figure.tight_layout(rect=[0, 0, 1, 0.96])

def save_data_file(filename, calculation_results, config):
    """Save detailed data file"""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"# Ion Radiation Calculation Results\n")
        f.write(f"# Generated: {datetime.now()}\n\n")
        
        f.write(f"# Ion: {config['ion']['name']} (Z={config['ion']['atomic_number']})\n")
        f.write(f"# Target Material: {config['target_material']['name']}\n")
        f.write(f"# Initial Energy: {config['ion']['energy_MeV_per_amu']} MeV/amu\n\n")
        
        f.write(f"# Input Parameters:\n")
        for section, params in config.items():
            if isinstance(params, dict):
                f.write(f"# {section}:\n")
                for key, value in params.items():
                    f.write(f"#   {key}: {value}\n")
        f.write(f"\n")
        
        if 'energy_info' in calculation_results:
            energy_info = calculation_results['energy_info']
            f.write(f"# SRIM Energy Loss Analysis:\n")
            for key, value in energy_info.items():
                f.write(f"# {key}: {value}\n")
            f.write(f"\n")
        
        analysis_info = calculation_results['analysis_info']
        f.write(f"# Analysis Results:\n")
        f.write(f"# Molecular density: {analysis_info['molecular_density_nm3']:.6e} molecules/nm³\n")
        f.write(f"# Track center temperature: {analysis_info['track_center_temperature_K']:,.0f} K\n")
        f.write(f"# Total deposited energy: {analysis_info['total_energy_keV']:.6e} keV\n\n")
        
        f.write("radius_nm\toriginal_dose_keV_nm3\tnormalized_dose_keV_nm3\ttemperature_K\tenergy_per_atom_eV\n")
        
        radius = calculation_results['radius']
        original_dose = calculation_results['original_dose']
        normalized_dose = calculation_results['normalized_dose']
        temperature = calculation_results['temperature']
        energy_per_atom = calculation_results['analysis_info']['energy_per_atom_eV']
        
        for r, d_o, d_n, t, e_a in zip(radius, original_dose, normalized_dose, temperature, energy_per_atom):
            f.write(f"{r:.6e}\t{d_o:.6e}\t{d_n:.6e}\t{t:.6e}\t{e_a:.6e}\n")

def load_config(config_path):
    """Load configuration from JSON file"""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file '{config_path}' not found.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in configuration file: {e}")
        sys.exit(1)