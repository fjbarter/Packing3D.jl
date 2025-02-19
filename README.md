# Packing3D.jl #
Packing3D.jl is a Julia package for calculating the local and bulk packing density of spherical particles in 3D. It supports both Cartesian and cylindrical coordinate systems and offers robust meshing capabilities for defining regions and visualising packing density distributions. A standout feature is the proprietary `read_vtk_file` function, which reliably parses legacy VTK files to extract particle coordinates and associated data. **Python version available [here](https://github.com/fjbarter/packing3d).**

![Example Image1](https://github.com/fjbarter/packing3d/blob/main/source/Before_After_Vibration.png?raw=true) *Cross-sectional packing density distribution of 100,000 particles (~500-600 microns) in a cylindrical container (75 mm diameter), computed over z = [0.005 m, 0.020 m].* 

## Table of Contents ##
- [Requirements](#requirements)
- [Installation](#installation)
- [Key Functions](#key-functions)
  - [`compute_packing_cartesian`](#compute_packing_cartesian)
  - [`compute_packing_cylindrical`](#compute_packing_cylindrical)
  - [`calculate_segregation_intensity`](#calculate_segregation_intensity)
  - [`calculate_lacey`](#calculate-lacey)
  - [`convert_to_cylindrical`](#convert_to_cylindrical)
  - [`read_vtk_file`](#read_vtk_file)
  - [`retrieve_coordinates`](#retrieve_coordinates)
  - [`get_mesh_bounds`](#get_mesh_bounds)
  - [`split_data` and `match_split_data`](#split_data-and-match_split_data)
  - [`extract_points`](#extract_points)
- [Mesh Structure](#mesh-structure)
- [How It Works](#how-it-works)
- [Testing](#testing)
- [Examples](#examples)
- [Limitations](#limitations)
- [Planned Features](#planned-features)
- [License](#license)
- [Contact](#contact)

--- 

## Requirements

- **Julia 1.6 or later**

Packing3D.jl is written entirely in Julia and uses standard libraries. No external dependencies are required beyond what is available in the Julia ecosystem.

--- 

## Installation

The package is not yet registered with Julia's package manager. To install, start Julia and run:
```
using Pkg
Pkg.develop(url="https://github.com/fjbarter/packing3d_jl")
```


--- 

## Key Functions

### `compute_packing_cartesian`

**Description:**  
Computes the packing density of particles within a user-defined Cartesian region. Particle data can be provided directly or loaded from a VTK file.

**Arguments:**
- `file::Union{String, Nothing}`: Path to the VTK file (if particle arrays are not provided).
- `x_data, y_data, z_data::Union{Vector{Float64}, Nothing}`: Particle coordinates.
- `radii::Union{Vector{Float64}, Nothing}`: Particle radii.
- `boundaries::Union{Vector{Float64}, Nothing}`: A 6-element array specifying `[x_min, x_max, y_min, y_max, z_min, z_max]`.

**Returns:**  
A `Float64` representing the fraction of the cell volume occupied by the particles.

**Example:**
using Packing3D

density = compute_packing_cartesian(file="particles.vtk")
println("Packing Density: ", density)

--- 

### `compute_packing_cylindrical`

**Description:**  
Calculates the packing density of particles within a cylindrical region, accounting for radial, angular, and axial boundaries.

**Arguments:**
- `file::Union{String, Nothing}`: Path to the VTK file.
- `r_data, theta_data, z_data::Union{Vector{Float64}, Nothing}`: Particle coordinates in cylindrical form.
- `radii::Union{Vector{Float64}, Nothing}`: Particle radii.
- `boundaries::Union{Vector{Float64}, Nothing}`: A 6-element array specifying `[r_min, r_max, theta_min, theta_max, z_min, z_max]`.
- `accurate_cylindrical::Bool`: Toggle for a more accurate (but computationally intensive) overlap calculation.

**Returns:**  
A `Float64` representing the packing density.

**Example:**
using Packing3D

density = compute_packing_cylindrical(file="particles.vtk")
println("Cylindrical Packing Density: ", density)

--- 

### `calculate_segregation_intensity`

**Description:**  
Evaluates the segregation intensity between two particle datasets (e.g., small vs. large particles) within a cylindrical mesh.

**Arguments:**
- `data_1::Dict`, `data_2::Dict`: Particle datasets.
- `cylinder_radius`, `cylinder_base_level`, `cylinder_height::Float64`: Cylindrical region parameters.
- `target_num_cells::Int`: Approximate number of cells for the mesh.
- Optional keyword arguments to control cell output and volume calculations.

**Returns:**  
A `Float64` representing the segregation intensity (0 for perfectly mixed, 1 for completely segregated).

--- 

### `calculate_lacey`

**Description:**  
Computes the Lacey mixing index between two particle datasets in a cylindrical mesh, providing a quantitative measure of mixing.

**Arguments:**  
Same as for `calculate_segregation_intensity`, with additional options for output and clamping.

**Returns:**  
A `Float64` representing the Lacey index (ranging from 0 for perfect mixing to 1 for complete segregation).

--- 

### `convert_to_cylindrical`

**Description:**  
Converts Cartesian coordinates to cylindrical coordinates with angles mapped to the range [0, 2π].

**Arguments:**
- `x_data, y_data::Vector{Float64}`

**Returns:**  
A tuple `(r_data, theta_data)`.

**Example:**
```
using Packing3D

r, theta = convert_to_cylindrical([1.0, 2.0], [1.0, 0.0])
println("r: ", r, " theta: ", theta)
```
--- 

### `read_vtk_file`

**Description:**  
A proprietary function that reads legacy VTK files (ASCII POLYDATA format) and extracts point coordinates along with any associated data (scalars or fields). This function is critical for robustly importing particle datasets into Packing3D.jl.

**Arguments:**
- `file::String`: Path to the `.vtk` file.

**Returns:**  
A `Dict` containing:
- `:points` – an N×3 matrix of point coordinates.
- `:point_data` – a dictionary of associated attributes.

**Example:**
```
using Packing3D

data = read_vtk_file("particles.vtk")
println("Number of points: ", size(data[:points], 1))
```
--- 

### `retrieve_coordinates`

**Description:**  
Extracts x, y, z coordinates and radii from the dataset produced by `read_vtk_file`.

**Returns:**  
A tuple `(x_data, y_data, z_data, radii)`.

--- 

### `get_mesh_bounds`

**Description:**  
Computes the axis-aligned bounding box for a mesh based on its point coordinates.

**Returns:**  
A vector `[x_min, x_max, y_min, y_max, z_min, z_max]`.

--- 

### `split_data` and `match_split_data`

**Description:**  
These functions allow you to partition your dataset into subsets (e.g., based on spatial criteria or particle properties).  
- `split_data` returns sets of point IDs for each subset.
- `match_split_data` extracts the corresponding subsets from a dataset.

--- 

### `extract_points`

**Description:**  
Extracts a subset of points and their associated data based on a boolean mask.

--- 

### Mesh Structure

The `Mesh` struct in Packing3D.jl represents a computational mesh in either Cartesian or cylindrical coordinates. It encapsulates:
- The coordinate system (`:cartesian` or `:cylindrical`).
- The number of divisions along each axis.
- The overall cell boundaries and indices.

**Example:**
```
using Packing3D

Example parameters for a cylindrical mesh:
divisions = Dict("r"=>5, "theta"=>3, "z"=>10)
params = Dict("cylinder_radius"=>10.0, "cylinder_base_level"=>0.0, "cylinder_height"=>20.0)
mesh = Mesh(:cylindrical, divisions; params=params)
println("Total mesh cells: ", get_total_cells(mesh))
```
--- 

## How It Works

1. **Data Input:**  
   The package relies on the proprietary `read_vtk_file` to import particle data from legacy VTK files. This function extracts both spatial coordinates and associated scalar or field data.

2. **Coordinate Conversion:**  
   Cartesian coordinates can be converted to cylindrical coordinates via `convert_to_cylindrical`, ensuring compatibility with cylindrical calculations.

3. **Boundary Determination & Meshing:**  
   Boundaries are either automatically computed (using provided particle data) or specified by the user. The `Mesh` struct (via the MeshModule) partitions the defined region into cells for localised density calculations.

4. **Overlap & Volume Calculation:**  
   For particles that partially intersect the defined boundaries, the package calculates the overlapping volumes using analytical and numerical integration techniques.

5. **Density & Mixing Metrics:**  
   Packing density is computed as the ratio of particle volume (including partial contributions) to cell volume. Additionally, segregation and mixing indices are provided to quantify the distribution of different particle types.

--- 

## Testing

- **Unit Tests:**  
  Comprehensive unit tests will be provided at a later date

--- 

## Examples

Sample scripts demonstrating usage are provided in the `examples` directory. These include:
- Calculating packing density in Cartesian and cylindrical systems.
- Splitting datasets based on spatial criteria.
- Computing segregation and Lacey mixing indices.

--- 

## Limitations

- **Particle Shape:**  
  The package is optimised for use with spherical particles.
  
- **VTK Format:**  
  Only legacy ASCII POLYDATA VTK files are supported with the inbuilt VTK reader. Files must be correctly formatted to be parsed by `read_vtk_file`.
  If you would like to parse XML VTK files, please see the **[ReadVTK.jl](https://github.com/JuliaVTK/ReadVTK.jl)** library.

--- 

## Planned Features

- **Visualisation:**  
  Integration with interactive plotting libraries for 3D visualisation of packing densities and meshes.
  
- **Extended Particle Types:**  
  Support for multisphere particles will be added in the near future.
  
- **Performance Optimisations:**  
  Easy parallel processing options would be nice to add in.

--- 

## License

This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

--- 

## Contact

For questions, support, or contributions, please contact:
- **Name:** Freddie Barter  
- **Email:** [fjbarter@outlook.com](mailto:fjbarter@outlook.com)
