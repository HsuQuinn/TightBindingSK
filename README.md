# SlaterKosterTB Project

**Author**: Hsu Quinn  
**Date**: March 2025  

## Project Overview
SlaterKosterTB is a Julia-based project for constructing and visualizing the electronic band structure and energy contour plots of crystalline materials using the Tight-Binding Model. The project leverages Slater-Koster coefficients to parameterize the tight-binding Hamiltonian, enabling accurate modeling of various crystal structures and orbital interactions. It is particularly suitable for studying electronic structures in solid-state physics.

## Features
- Define crystal structures and atomic orbitals.
- Generate Brillouin zone paths and compute band structures.
- Plot band structure diagrams and energy contour plots.
- Support multiple orbital types (e.g., `s`, `p`, `d` orbitals).

## Directory Structure
- `src/`: Source code of the project, including the following modules:
  - `structure.jl`: Defines crystal structures and atoms.
  - `klib.jl`: Generates Brillouin zone paths.
  - `hamilton.jl`: Generates Hamiltonians and computes band structures.
  - `param.jl`: Defines parameters for the tight-binding model.
  - `draw.jl`: Plots band structures and energy contour plots.
- `test/`: Test codes with examples of different crystal structures (e.g., graphene, honeycomb structures).
- `README.md`: Project documentation.
- `.gitignore`: Git ignore file.

## Installation
Run the following commands to install the required Julia packages:
```julia
using Pkg
Pkg.add("Plots")
Pkg.add("LinearAlgebra")
```

## Usage
1. Clone the repository:
   ```bash
   git clone <repository-url>
   ```
2. Navigate to the project directory:
   ```bash
   cd SKTB
   ```
3. Run the test scripts:
   ```bash
   julia test/test_graphene.jl
   ```
   Or run other test files (e.g., `test/test_honey.jl`).

## Example
Here is an example of running `test/test_graphene.jl`:
1. Define the crystal structure and orbitals for graphene.
2. Compute the band structure along the Brillouin zone path.
3. Plot the band structure and energy contour diagrams, which are saved as PNG files.

## Output
After running the test scripts, the generated image files will be saved in the current directory, such as:
- `band_graphene.png`: Band structure diagram for graphene.
- `contour_graphene.png`: Energy contour plot for graphene.

## Contribution
We welcome issues and contributions! Please follow these steps:
1. Fork the repository.
2. Create a new branch and commit your changes.
3. Submit a Pull Request.

## License
This project is licensed under the MIT License.

## TODO List
Planned updates for future versions:
1. Implement orbital projection.
2. Add Spin-Orbit Coupling (SOC).
3. Include next-nearest-neighbor interactions.
4. Simulate ARPES (Angle-Resolved Photoemission Spectroscopy).
