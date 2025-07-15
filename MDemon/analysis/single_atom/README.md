# Single Atom RDF Analysis

This module provides radial distribution function (RDF) analysis capabilities for single atoms in MDemon.

## Features

- **Basic RDF Analysis**: Calculate RDF for individual atoms or groups of atoms
- **Time-Resolved RDF**: Analyze RDF evolution over trajectory frames
- **Parallel Processing**: Built-in Dask-based parallel processing for large systems
- **Atom Selection**: Flexible atom selection using boolean masks
- **Visualization**: Built-in plotting capabilities

## Quick Start

### Basic RDF Analysis

```python
from MDemon.analysis.single_atom import RDFAnalyzer
from MDemon.core import Universe

# Load your system
universe = Universe("trajectory.dump", format="LAMMPS")

# Create RDF analyzer
analyzer = RDFAnalyzer(
    universe=universe,
    r_range=(0.0, 10.0),  # Distance range in Angstroms
    n_bins=100,           # Number of bins
    scheduler='threads'   # Parallel scheduler
)

# Analyze single atom
atom_index = 0
r_values, rdf_values = analyzer.analyze_single_atom(atom_index)

# Analyze multiple atoms in parallel
atom_indices = [0, 1, 2, 3, 4]
result = analyzer.analyze_parallel_rdf(atom_indices=atom_indices)

# Get results for specific atom
r_vals, rdf_vals = result.get_rdf(atom_index=0)

# Plot results
result.plot()
```

### Atom Selection

```python
import numpy as np

# Select atoms by species
species_mask = np.array([atom.species == 1 for atom in universe.atoms])

# Select atoms by spatial region
z_mask = np.array([atom.coordinate[2] > 5.0 for atom in universe.atoms])

# Combine selections
combined_mask = species_mask & z_mask

# Create analyzer with selection
analyzer = RDFAnalyzer(
    universe=universe,
    atom_selection=combined_mask,
    r_range=(0.0, 10.0),
    n_bins=100
)

# Analyze selected atoms
result = analyzer.analyze_parallel_rdf()
```

### Time-Resolved RDF

```python
from MDemon.analysis.single_atom import TimeResolvedRDFAnalyzer

# Create time-resolved analyzer
analyzer = TimeResolvedRDFAnalyzer(
    universe=universe,
    r_range=(0.0, 10.0),
    n_bins=100,
    time_window=(0, 100),  # Frame range
    scheduler='threads'
)

# Analyze trajectory
atom_indices = [0, 1, 2]
trajectory_results = analyzer.analyze_trajectory(
    atom_indices=atom_indices,
    frame_step=2  # Every 2nd frame
)

# Compute time-averaged RDF
r_values, time_avg_rdf = analyzer.compute_time_averaged_rdf(atom_indices)

# Compute time correlation
correlation = analyzer.compute_time_correlation(trajectory_results)
```

## Performance Tips

1. **Use appropriate scheduler**:
   - `'threads'`: Good for I/O-bound tasks
   - `'processes'`: Better for CPU-intensive calculations
   - `'distributed'`: For cluster computing

2. **Optimize bin count**: More bins = higher resolution but slower computation

3. **Use atom selection**: Analyze only relevant atoms to save time

4. **Batch processing**: For large systems, the analyzer automatically optimizes memory usage

## Running Tests

```bash
# Run all RDF tests
uv run python -m pytest tests/test_single_atom_rdf.py -v

# Run specific test
uv run python -m pytest tests/test_single_atom_rdf.py::TestRDFAnalyzer::test_initialization -v
```

## Examples

See `examples/rdf_example.py` for comprehensive usage examples including:
- Basic RDF analysis
- Time-resolved analysis
- Atom selection techniques
- Visualization

Run the examples:
```bash
uv run python examples/rdf_example.py
```

## API Reference

### RDFAnalyzer

Main class for RDF analysis.

**Parameters:**
- `universe`: MDemon Universe object
- `atom_selection`: Boolean mask for atom selection (optional)
- `r_range`: Distance range tuple (r_min, r_max)
- `n_bins`: Number of distance bins
- `scheduler`: Dask scheduler ('threads', 'processes', 'distributed')

**Methods:**
- `analyze_single_atom(atom_index)`: Analyze single atom
- `analyze_parallel_rdf(atom_indices)`: Parallel analysis of multiple atoms
- `compute_average_rdf(atom_indices)`: Compute average RDF

### TimeResolvedRDFAnalyzer

Extends RDFAnalyzer for time-dependent analysis.

**Additional Parameters:**
- `time_window`: Frame range tuple (start, end)

**Methods:**
- `analyze_trajectory(atom_indices, frame_step)`: Analyze trajectory
- `compute_time_averaged_rdf(atom_indices)`: Time-averaged RDF
- `compute_time_correlation(rdf_series)`: Time correlation function

### RDFResult

Result container for RDF analysis.

**Methods:**
- `get_rdf(atom_index)`: Get RDF for specific atom
- `get_first_peak(atom_index)`: Get first peak position
- `plot(atom_index)`: Plot RDF

## Dependencies

- numpy >= 1.24.0
- dask[complete] >= 2023.0.0
- matplotlib >= 3.6.0 (for plotting)
- MDemon core modules

## License

This module is part of MDemon and follows the same license terms. 