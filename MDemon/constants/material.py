"""
Material Properties Constants

Common material properties used in molecular dynamics simulations.
Particularly useful for materials related to Ga₂O₃ research and general MD simulations.
"""


class MaterialConstants:
    """Container for material properties constants"""
    
    # Common material densities [kg/m³]
    densities = {
        # Metals
        'Al': 2700.0,
        'Cu': 8960.0,
        'Fe': 7874.0,
        'Au': 19300.0,
        'Ag': 10490.0,
        'Pt': 21450.0,
        'Ti': 4506.0,
        'Ni': 8908.0,
        
        # Semiconductors
        'Si': 2329.0,
        'Ge': 5323.0,
        'GaAs': 5316.0,
        'InP': 4810.0,
        'GaN': 6150.0,
        
        # Oxides (relevant for Ga₂O₃ research)
        'Ga2O3_beta': 5950.0,  # β-Ga₂O₃
        'Ga2O3_alpha': 6440.0,  # α-Ga₂O₃
        'Al2O3': 3980.0,  # Sapphire
        'SiO2': 2650.0,  # Quartz
        'TiO2': 4230.0,  # Rutile
        'ZnO': 5606.0,
        'MgO': 3580.0,
        'CaO': 3340.0,
        
        # Common liquids
        'water': 1000.0,
        'ethanol': 789.0,
        'benzene': 876.0,
        'acetone': 784.0,
        
        # Polymers
        'polyethylene': 920.0,
        'polystyrene': 1050.0,
        'PVC': 1380.0,
        'PMMA': 1180.0,
    }
    
    # Melting points [K]
    melting_points = {
        # Metals
        'Al': 933.47,
        'Cu': 1357.77,
        'Fe': 1811.0,
        'Au': 1337.33,
        'Ag': 1234.93,
        'Pt': 2041.4,
        'Ti': 1941.0,
        'Ni': 1728.0,
        
        # Semiconductors
        'Si': 1687.0,
        'Ge': 1211.4,
        'GaN': 2790.0,  # Decomposes
        
        # Oxides
        'Ga2O3_beta': 2123.0,  # β-Ga₂O₃ (decomposes)
        'Al2O3': 2345.0,
        'SiO2': 1996.0,
        'TiO2': 2116.0,
        'ZnO': 2248.0,  # Decomposes
        'MgO': 3125.0,
        
        # Common substances
        'water': 273.15,
        'ethanol': 159.05,
        'benzene': 278.68,
    }
    
    # Thermal conductivity [W/(m⋅K)]
    thermal_conductivity = {
        # Metals (room temperature)
        'Al': 237.0,
        'Cu': 401.0,
        'Fe': 80.4,
        'Au': 318.0,
        'Ag': 429.0,
        'Pt': 71.6,
        'Ti': 21.9,
        'Ni': 90.9,
        
        # Semiconductors
        'Si': 149.0,
        'Ge': 60.2,
        'GaN': 130.0,
        
        # Oxides
        'Ga2O3_beta': 27.0,  # β-Ga₂O₃ (along b-axis)
        'Al2O3': 35.0,  # Sapphire (perpendicular to c-axis)
        'SiO2': 1.38,  # Quartz
        'TiO2': 11.8,  # Rutile
        'ZnO': 54.0,
        'MgO': 42.0,
        
        # Liquids
        'water': 0.606,
        'ethanol': 0.167,
        'benzene': 0.141,
    }
    
    # Bulk modulus [GPa]
    bulk_modulus = {
        # Metals
        'Al': 76.0,
        'Cu': 140.0,
        'Fe': 170.0,
        'Au': 180.0,
        'Ag': 100.0,
        'Pt': 230.0,
        'Ti': 110.0,
        'Ni': 180.0,
        
        # Semiconductors
        'Si': 97.6,
        'Ge': 75.8,
        'GaN': 210.0,
        
        # Oxides
        'Ga2O3_beta': 220.0,  # β-Ga₂O₃ (estimated)
        'Al2O3': 252.0,
        'SiO2': 37.0,
        'TiO2': 211.0,
        'ZnO': 142.6,
        'MgO': 162.0,
    }
    
    # Young's modulus [GPa]
    youngs_modulus = {
        # Metals
        'Al': 70.0,
        'Cu': 110.0,
        'Fe': 211.0,
        'Au': 78.0,
        'Ag': 83.0,
        'Pt': 168.0,
        'Ti': 116.0,
        'Ni': 200.0,
        
        # Semiconductors
        'Si': 112.4,  # <100> direction
        'Ge': 102.7,
        'GaN': 295.0,
        
        # Oxides
        'Ga2O3_beta': 230.0,  # β-Ga₂O₃ (estimated)
        'Al2O3': 400.0,
        'SiO2': 73.0,
        'TiO2': 283.0,
        'ZnO': 111.2,
        'MgO': 248.0,
    }
    
    # Lattice parameters for common crystal structures [Å]
    lattice_parameters = {
        # Cubic metals (fcc, bcc)
        'Al_fcc': {'a': 4.046},
        'Cu_fcc': {'a': 3.615},
        'Fe_bcc': {'a': 2.867},
        'Au_fcc': {'a': 4.078},
        'Ag_fcc': {'a': 4.085},
        'Pt_fcc': {'a': 3.924},
        'Ni_fcc': {'a': 3.524},
        
        # Diamond/zincblende structures
        'Si_diamond': {'a': 5.431},
        'Ge_diamond': {'a': 5.658},
        'GaAs_zb': {'a': 5.653},
        
        # Hexagonal structures
        'GaN_wurtzite': {'a': 3.189, 'c': 5.185},
        'ZnO_wurtzite': {'a': 3.250, 'c': 5.207},
        
        # Corundum structure (α-Al₂O₃, α-Ga₂O₃)
        'Al2O3_corundum': {'a': 4.759, 'c': 12.991},
        'Ga2O3_alpha': {'a': 4.982, 'c': 13.433},
        
        # β-Ga₂O₃ (monoclinic)
        'Ga2O3_beta': {'a': 12.214, 'b': 3.037, 'c': 5.798, 'beta': 103.83},
        
        # Rocksalt structures
        'MgO_rocksalt': {'a': 4.211},
        'CaO_rocksalt': {'a': 4.815},
    }
    
    @classmethod
    def get_density(cls, material):
        """Get density for a material [kg/m³]"""
        return cls.densities.get(material, None)
    
    @classmethod
    def get_melting_point(cls, material):
        """Get melting point for a material [K]"""
        return cls.melting_points.get(material, None)
    
    @classmethod
    def get_thermal_conductivity(cls, material):
        """Get thermal conductivity for a material [W/(m⋅K)]"""
        return cls.thermal_conductivity.get(material, None)
    
    @classmethod
    def get_bulk_modulus(cls, material):
        """Get bulk modulus for a material [GPa]"""
        return cls.bulk_modulus.get(material, None)
    
    @classmethod
    def get_youngs_modulus(cls, material):
        """Get Young's modulus for a material [GPa]"""
        return cls.youngs_modulus.get(material, None)
    
    @classmethod
    def get_lattice_parameters(cls, structure):
        """Get lattice parameters for a crystal structure [Å]"""
        return cls.lattice_parameters.get(structure, None)
    
    @classmethod
    def list_available_materials(cls):
        """List all available materials"""
        materials = set()
        for prop_dict in [cls.densities, cls.melting_points, 
                         cls.thermal_conductivity, cls.bulk_modulus, 
                         cls.youngs_modulus]:
            materials.update(prop_dict.keys())
        return sorted(list(materials))


# Create the singleton instance
MATERIAL = MaterialConstants() 