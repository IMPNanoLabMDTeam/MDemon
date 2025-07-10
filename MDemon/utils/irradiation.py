import json
import math
import sys
from datetime import datetime

import numpy as np

from ..constants import CONVERSION, IRRADIATION, MATERIAL, PHYSICS


class WaligorskiZhangCalculator:
    """
    Waligorski-Zhang radial dose distribution model calculator

    This calculator uses constants from the MDemon.constants module and
    requires material properties to be specified for calculations.
    """

    def __init__(
        self,
        molecular_weight,
        atoms_per_molecule,
        density_g_per_cm3,
        ion_name,
        ion_energy_MeV_per_amu,
        ion_Z,
        energy_loss_keV_per_um,
        g_factor=0.17,
    ):
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

        # irradiation parameters
        self.g_factor = g_factor

    def calculate_beta_gamma(self, energy_MeV_per_amu):
        """Calculate relativistic β and γ factors"""
        γ = 1 + energy_MeV_per_amu / self.amu_to_MeV
        β = math.sqrt(1 - 1 / γ**2)
        return β, γ

    def calculate_effective_charge(self, Z, β):
        """Calculate effective charge using Barkas correction"""
        return Z * (1 - math.exp(-self.barkas_coefficient * β * Z ** (-2 / 3)))

    def calculate_range_parameter_alpha(self, β):
        """Calculate range parameter α"""
        return IRRADIATION.get_range_parameter_alpha(β)

    def calculate_radial_dose(self, radius_nm):
        """Calculate radial dose distribution"""
        β, γ = self.calculate_beta_gamma(self.ion_energy_MeV_per_amu)
        Z_eff = self.calculate_effective_charge(self.ion_Z, β)
        α = self.calculate_range_parameter_alpha(β)

        ionization_energy_keV = self.ionization_energy_eV / 1000.0
        θ = self.k_g_per_cm2_keV_per_alpha * (ionization_energy_keV**1.079)
        W = 2 * self.electron_mass_keV * β**2 * γ**2
        T = self.k_g_per_cm2_keV_per_alpha * (W**α)

        radius_cm = radius_nm * 1e-7
        t_g_per_cm2 = radius_cm * self.density_g_per_cm3

        if t_g_per_cm2 <= 0 or (T + θ) <= 0 or α <= 0:
            return 0

        water_constant_keV_per_cm = self.water_constant_keV_per_mm * 1e1
        constant_factor = (
            water_constant_keV_per_cm * Z_eff**2 / (2 * np.pi * α * β**2 * t_g_per_cm2)
        )
        fraction = (t_g_per_cm2 + θ) / (T + θ)

        if fraction >= 1:
            return 0

        power_term = (1 - fraction) ** (1 / α)
        dose_keV_cm3_per_g2 = constant_factor * power_term / (t_g_per_cm2 + θ)
        dose_keV_per_cm3 = dose_keV_cm3_per_g2 * (self.density_g_per_cm3) ** 2
        dose_keV_per_nm3 = dose_keV_per_cm3 * 1e-21
        # 一次性输出所有变量用于debug
        # print(f"dose_keV_per_nm3: {dose_keV_per_nm3: .3e}, constant_factor: {constant_factor: .3e}, power_term: {power_term: .3e}, t_g_per_cm2: {t_g_per_cm2: .3e}, density_g_per_cm3: {self.density_g_per_cm3: .3e}, ionization_energy_keV: {ionization_energy_keV: .3e}, θ: {θ: .3e}, T: {T: .3e}, α: {α: .3e}, β: {β: .3e}, γ: {γ: .3e}, Z_eff: {Z_eff: .3e}, W: {W: .3e}")
        return dose_keV_per_nm3

    def calculate_radial_energy_distribution(self, radius_array_nm):
        """Calculate complete radial energy distribution"""
        dose_array = np.array([self.calculate_radial_dose(r) for r in radius_array_nm])

        atomic_density_nm3 = (
            (self.density_g_per_cm3 * self.avogadro_constant * 1e-21)
            / self.molecular_weight
            * self.atoms_per_molecule
        )

        # Calculate cumulative energy
        cumulative_energy = np.zeros_like(radius_array_nm)
        for i in range(1, len(radius_array_nm)):
            dr = radius_array_nm[i] - radius_array_nm[i - 1]
            ring_area = 2 * np.pi * radius_array_nm[i] * dr
            ring_energy = dose_array[i] * ring_area * 1.0
            cumulative_energy[i] = cumulative_energy[i - 1] + ring_energy

        # normalize the dose by the energy loss per um
        normalized_dose = (
            dose_array
            * self.g_factor
            * self.energy_loss_keV_per_um
            / cumulative_energy[-1]
            / 1000
        )

        # calculate the temperature
        energy_per_atom_eV = (normalized_dose * 1000.0) / atomic_density_nm3
        temperature = energy_per_atom_eV * self.eV_to_K

        self.results = {
            "radius": radius_array_nm,
            "dose_density": dose_array,
            "normalized_dose_density": normalized_dose,
            "energy_per_atom": energy_per_atom_eV,
            "temperature": temperature,
            "cumulative_energy": cumulative_energy,
            "atomic_density": atomic_density_nm3,
            "calculated_energy_loss": (
                cumulative_energy[-1] if len(cumulative_energy) > 0 else 0
            ),
        }
        # print(f"dose_density: {self.results['dose_density']}")
        return self.results

    def save_results(self, filepath=None, include_metadata=True):
        """
        Save calculation results to a JSON file

        Parameters:
        -----------
        filepath : str, optional
            Path to save the JSON file. If None, generates a default filename based on
            ion properties and timestamp. Default is None.
        include_metadata : bool, optional
            Whether to include calculation parameters as metadata. Default is True.

        Returns:
        --------
        bool
            True if save was successful, False otherwise
        """
        if not hasattr(self, "results") or self.results is None:
            print(
                "No results to save. Please run calculate_radial_energy_distribution() first."
            )
            return False

        # Generate default filepath if not provided
        if filepath is None:
            import os

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = (
                f"{self.ion_name}_{self.ion_energy_MeV_per_amu}MeVamu_{timestamp}.json"
            )

            # Create results directory if it doesn't exist
            results_dir = "results"
            if not os.path.exists(results_dir):
                os.makedirs(results_dir)

            filepath = os.path.join(results_dir, filename)

        try:
            # Prepare data for JSON serialization
            save_data = {}

            # Add metadata if requested
            if include_metadata:
                save_data["metadata"] = {
                    "ion_name": self.ion_name,
                    "ion_energy_MeV_per_amu": self.ion_energy_MeV_per_amu,
                    "ion_Z": self.ion_Z,
                    "energy_loss_keV_per_um": self.energy_loss_keV_per_um,
                    "g_factor": self.g_factor,
                    "molecular_weight": self.molecular_weight,
                    "atoms_per_molecule": self.atoms_per_molecule,
                    "density_g_per_cm3": self.density_g_per_cm3,
                    "timestamp": datetime.now().isoformat(),
                    "calculator": "WaligorskiZhangCalculator",
                }

            # Add results data (convert numpy arrays to lists for JSON serialization)
            save_data["results"] = {}
            for key, value in self.results.items():
                if isinstance(value, np.ndarray):
                    save_data["results"][key] = value.tolist()
                else:
                    save_data["results"][key] = value

            # Save to JSON file
            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(save_data, f, indent=2, ensure_ascii=False)

            print(f"Results successfully saved to: {filepath}")
            return True

        except Exception as e:
            print(f"Error saving results to {filepath}: {e}")
            return False

    def plot_results(
        self, save_path=None, show_plot=True, log_scale=True, xlim=None, ylim=None
    ):
        """
        Plot calculation results

        Parameters:
        -----------
        save_path : str, optional
            Path to save the plot. If None, plot is not saved.
        show_plot : bool, optional
            Whether to display the plot. Default is True.
        log_scale : bool, optional
            Whether to use logarithmic scale for plots. If True, uses log-log or semi-log scales.
            If False, uses linear scales. Default is True.
        xlim : tuple, optional
            X-axis limits as (xmin, xmax). Applied to all subplots. If None, uses automatic scaling.
        ylim : dict, optional
            Y-axis limits for each subplot. Keys: 'dose', 'energy', 'temperature', 'cumulative'.
            Values: tuples (ymin, ymax). If None or key missing, uses automatic scaling.
            Example: {'dose': (1e-10, 1e-5), 'temperature': (100, 10000)}

        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print(
                "Warning: matplotlib is required for plotting. Please install it with 'pip install matplotlib'"
            )
            return None

        if not hasattr(self, "results") or self.results is None:
            print(
                "No results to plot. Please run calculate_radial_energy_distribution() first."
            )
            return None

        # Create figure with subplots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(
            f"Waligorski-Zhang Model Results for {self.ion_name} ion\n"
            f"Energy: {self.ion_energy_MeV_per_amu} MeV/amu, dE/dx: {self.energy_loss_keV_per_um} keV/μm",
            fontsize=14,
            fontweight="bold",
        )

        radius = self.results["radius"]

        # Plot 1: Dose density distribution
        if log_scale:
            ax1.loglog(
                radius,
                self.results["dose_density"],
                "b-",
                linewidth=2,
                label="Original dose",
            )
            ax1.loglog(
                radius,
                self.results["normalized_dose_density"],
                "r--",
                linewidth=2,
                label="Normalized dose",
            )
        else:
            ax1.plot(
                radius,
                self.results["dose_density"],
                "b-",
                linewidth=2,
                label="Original dose",
            )
            ax1.plot(
                radius,
                self.results["normalized_dose_density"],
                "r--",
                linewidth=2,
                label="Normalized dose",
            )
        ax1.set_xlabel("Radius (nm)")
        ax1.set_ylabel("Dose Density (keV/nm³)")
        ax1.set_title("Radial Dose Distribution")
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        # Set axis limits
        if xlim is not None:
            ax1.set_xlim(xlim)
        if ylim is not None and "dose" in ylim:
            ax1.set_ylim(ylim["dose"])

        # Plot 2: Energy per atom
        if log_scale:
            ax2.loglog(radius, self.results["energy_per_atom"], "g-", linewidth=2)
        else:
            ax2.plot(radius, self.results["energy_per_atom"], "g-", linewidth=2)
        ax2.set_xlabel("Radius (nm)")
        ax2.set_ylabel("Energy per Atom (eV)")
        ax2.set_title("Energy Deposition per Atom")
        ax2.grid(True, alpha=0.3)
        # Set axis limits
        if xlim is not None:
            ax2.set_xlim(xlim)
        if ylim is not None and "energy" in ylim:
            ax2.set_ylim(ylim["energy"])

        # Plot 3: Temperature distribution
        if log_scale:
            ax3.loglog(radius, self.results["temperature"], "orange", linewidth=2)
        else:
            ax3.plot(radius, self.results["temperature"], "orange", linewidth=2)
        ax3.set_xlabel("Radius (nm)")
        ax3.set_ylabel("Temperature (K)")
        ax3.set_title("Temperature Distribution")
        ax3.grid(True, alpha=0.3)
        # Set axis limits
        if xlim is not None:
            ax3.set_xlim(xlim)
        if ylim is not None and "temperature" in ylim:
            ax3.set_ylim(ylim["temperature"])

        # Plot 4: Cumulative energy
        if log_scale:
            ax4.semilogx(radius, self.results["cumulative_energy"], "m-", linewidth=2)
        else:
            ax4.plot(radius, self.results["cumulative_energy"], "m-", linewidth=2)
        ax4.set_xlabel("Radius (nm)")
        ax4.set_ylabel("Cumulative Energy (keV·nm)")
        ax4.set_title("Cumulative Energy vs Radius")
        if xlim is not None:
            ax4.set_xlim(xlim)
        if ylim is not None and "cumulative" in ylim:
            ax4.set_ylim(ylim["cumulative"])
        ax4.grid(True, alpha=0.3)

        # Add information text box
        info_text = (
            f"Material: MW={self.molecular_weight} g/mol, "
            f"ρ={self.density_g_per_cm3} g/cm³\n"
            f"Ion: Z={self.ion_Z}, g-factor={self.g_factor}\n"
            f"Atomic density: {self.results['atomic_density']:.2e} nm⁻³\n"
            f"Calculated energy loss: {self.results['calculated_energy_loss']:.2f} keV·nm"
        )

        fig.text(
            0.02,
            0.02,
            info_text,
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
        )

        plt.tight_layout()

        # Save plot if path is provided
        if save_path:
            try:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                print(f"Plot saved to: {save_path}")
            except Exception as e:
                print(f"Warning: Could not save plot to {save_path}. Error: {e}")

        # Show plot if requested
        if show_plot:
            plt.show()

        return fig
