# Packing3D.jl #

[![Julia](https://github.com/fjbarter/Packing3D.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/fjbarter/Packing3D.jl/actions/workflows/ci.yml)

Jupyter Notebook available with Binder!

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fjbarter/Packing3D.jl/main?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Fpacking_and_mixing_demo.ipynb)

Packing3D.jl is a Julia package for calculating the local and bulk packing density of spherical particles in 3D. It supports both Cartesian and cylindrical coordinate systems and provides robust meshing capabilities for defining computational regions, along with tools for quantifying mixing and segregation of particle types. **Python version available [here](https://github.com/fjbarter/packing3d).**

![Example Image1](https://github.com/fjbarter/packing3d/blob/main/source/Before_After_Vibration.png?raw=true) *Cross-sectional packing density distribution of 100,000 particles (~500-600 microns) in a cylindrical container (75 mm diameter), computed over z = [0.005 m, 0.020 m].* 

**Public API Overview:**  
The package exposes a clean, high-level interface designed for users (such as Master's students) who may not be familiar with the underlying Julia intricacies. The main functions automatically dispatch to the appropriate Cartesian or cylindrical routines based on the provided `system` parameter.

The key public functions are:
- [`calculate_packing`](#calculate_packing)
- [`calculate_segregation_intensity`](#calculate_segregation_intensity)
- [`calculate_lacey`](#calculate_lacey)
- [`compute_volumes_per_cell`](#compute_volumes_per_cell)
- [`get_cell_index`](#get_cell_index)

In addition, the following functions and types are exported for public use:
- **Mesh Utilities:** `Mesh`, `get_mesh_bounds`, `get_total_cells`, `get_cell_volume`
- **Coordinate & I/O Functions:** `convert_to_cylindrical`, `read_vtk_file`, `save_vtk_file`, `retrieve_coordinates`
- **Data Partitioning Tools:** `split_data`, `match_split_data`, `extract_points`

> **Note:** Some internal docstrings may be outdated. The documentation below reflects the current intended usage.

## Table of Contents
- [Public API Overview](#public-api-overview)
- [Key Functions](#key-functions)
  - [`calculate_packing`](#calculate_packing)
  - [`calculate_segregation_intensity`](#calculate_segregation_intensity)
  - [`calculate_lacey`](#calculate_lacey)
  - [`compute_volumes_per_cell`](#compute_volumes_per_cell)
  - [`get_cell_index`](#get_cell_index)
  - [Mesh Utilities](#mesh-utilities)
  - [Coordinate Conversion & I/O](#coordinate-conversion--io)
  - [Data Partitioning Tools](#data-partitioning-tools)
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Examples](#examples)
- [Limitations](#limitations)
- [Planned Features](#planned-features)
- [Licence](#license)
- [Contact](#contact)

## Public API Overview

Packing3D.jl is built around a small set of high-level functions that abstract away the details of coordinate systems and meshing. Whether your data is in Cartesian or cylindrical form, the package automatically routes calls to the correct internal routines.

## Key Functions



Below are descriptions for the key public functions. In these descriptions, positional arguments and keyword arguments are distinguished where appropriate. Note that for many functions the mesh can be provided explicitly; if omitted, a mesh will be generated automatically based on the provided parameters.

---

### `calculate_packing`

**Description:**  
Calculates the packing density of particles within a defined region. This function supports both Cartesian and cylindrical systems and can load particle data from a VTK file if a direct data dictionary is not provided. Full cylindrical cells can be specified by either defining a negative `:r_min` or by using full angular boundaries (e.g. `[0, 2*pi]`). If a mesh is not supplied, one is generated automatically based on parameters such as the target number of cells.

**Positional and Keyword Arguments:**

- **Data Input:**  
  - `file::Union{String, Nothing}` (keyword): Path to a VTK file containing particle data, used if `data` is not provided.  
  - `data::Union{Dict, Nothing}` (keyword): A dictionary with particle data (if already loaded).

- **Boundaries:**  
  - For Cartesian systems:  
    - `boundaries`: A 6-element array `[x_min, x_max, y_min, y_max, z_min, z_max]` **or** a dictionary with keys `:x_min, :x_max, :y_min, :y_max, :z_min, :z_max`.  
  - For Cylindrical systems:  
    - `boundaries`: A 6-element array `[r_min, r_max, theta_min, theta_max, z_min, z_max]` **or** a dictionary with keys `:r_min, :r_max, :theta_min, :theta_max, :z_min, :z_max`.

- **System Specification:**  
  - `system::Symbol` (keyword): Specifies `:cartesian` or `:cylindrical`.

- **Additional Options:**  
  - Other system-specific parameters (e.g. `cylinder_radius`, `cylinder_base_level`, `cylinder_height`) may be supplied via keyword arguments.  
  - Any additional options (such as toggles for accurate overlap calculation) are also passed as keywords.

**Returns:**  
A `Float64` representing the fraction of the region's volume occupied by the particles.

**Example 1:** (Cylindrical)
```
using Packing3D

boundaries = Dict(
  :r_min => -1,  # Negative value indicates a full cylindrical cell
  :r_max => 0.0275,
  :theta_min => 0.0,
  :theta_max => 2*pi,
  :z_min => 0.01,
  :z_max => 0.04
)

packing_density = calculate_packing(
  ;
  file = "particles.vtk",
  boundaries = boundaries,
  system = :cylindrical
)

println("Packing Density: ", packing_density)
```

**Example 2:** (Cartesian)
```
using Packing3D

particles_file = "particles.vtk"
data = read_vtk_file(particles_file)
x_data, y_data, z_data, radii = retrieve_coordinates(data)

boundaries = Dict(
  :x_min => -0.03,
  :x_max => 0.03,
  :y_min => -0.03,
  :y_max => 0.03,
  :z_min => 0.01,
  :z_max => 0.04
)

packing_density = calculate_packing(
  ;
  x_data = x_data,
  y_data = y_data,
  z_data = z_data,
  radii = radii,
  boundaries = boundaries,
  system = :cartesian
)

println("Packing Density: ", packing_density)
```
---

### `calculate_lacey`

**Description:**  
Computes the Lacey mixing index between two particle datasets, providing a quantitative measure of mixing. It selects the appropriate internal routine based on the specified `system` and supports optional mesh input. Additional options allow for controlling output verbosity and for clamping the resulting index.

**Positional and Keyword Arguments:**

- **Data Input:**  
  - `data_1::Dict`, `data_2::Dict`: (Positional) Particle datasets for the two groups.

- **Mesh Input:**  
  - `mesh::Union{Mesh, Nothing}` (keyword): Optional mesh object. If not provided, a mesh is generated using `target_num_cells` and other parameters.

- **System & Mesh Parameters:**  
  - `system::Symbol` (keyword): Specify `:cartesian` or `:cylindrical`.  
  - `params::Dict` (keyword): for cylindrical systems, required parameters are `:cylinder_radius, :cylinder_base_level, :cylinder_height`. For Cartesian, `:x_min, :x_max, :y_min, :y_max, :z_min, :z_max`. Note that dictionary keys must be Symbols, denoted by a colon `:`.

- **Mesh Resolution:**  
  - `target_num_cells::Real` (keyword): Approximate number of mesh cells.

- **Optional Flags:**  
  - Flags to enable partial volume calculations (`calculate_partial_volumes`), verbosity (`verbose`), clamping of the final index (`clamp_0_to_1`), and cell number output (`output_num_cells`).  

**Returns:**  
A `Float64` representing the Lacey mixing index, where 1 indicates perfect mixing and 0 indicates complete segregation.

**Example 1**   
Providing data and `params` for a cylindrical mesh.
```
using Packing3D

params = Dict(
  :cylinder_radius => 0.03,
  :cylinder_base_level => 0.0,
  :cylinder_height => 0.08
)

lacey_index = calculate_lacey(
  data_1, data_2;
  system=:cylindrical,
  params=params,
  target_num_cells = 1000,
  calculate_partial_volumes = true,
  verbose = false,
  clamp_0_to_1 = false
)

println("Lacey Mixing Index: ", lacey_index)
```

**Example 2:**  
Providing a custom mesh with specified divisions in Cartesian coordinates (behaviour is identical in cylindrical).

```
using Packing3D

params = Dict(
  :x_min => -0.03,
  :x_max => 0.03,
  :y_min => -0.03,
  :y_max => 0.03,
  :z_min => 0.0,
  :z_max => 0.8
)

divisions => Dict(
  :x => 10,
  :y => 10,
  :z => 8
)

mesh = Mesh(:cartesian; params=params, divisions=divisions)

lacey_index = calculate_lacey(
  data_1, data_2;
  system = :cartesian,
  mesh = mesh
)

println("Lacey Mixing Index: ", lacey_index)
```

---

### `calculate_segregation_intensity`

**Description:**  
Functionality and arguments are identical to [`calculate_lacey`](#calculate_lacey)

**Returns:**  
A `Float64` representing the segregation intensity (0 for perfectly mixed, 1 for completely segregated).

---

### `compute_volumes_per_cell`

**Description:**  
Computes the volume contributions of particles to each cell in the computational mesh. This includes handling particles that partially overlap cell boundaries by calculating their partial volumes. This function is key for accurately determining local packing density. The function is the volume calculation kernel for `calculate_lacey` and `calculate_segregation_intensity`, but it is useful to have public access for custom analysis (e.g. local packing density calculation).

**Positional and Keyword Arguments:**

- **Data Input:**  
  - `data_1::Dict`, `data_2::Dict`: (Positional) Dictionaries containing particle data for each group.

- **Mesh Input:**  
  - `mesh::Mesh` (keyword): A mesh object defining the computational grid. This must be provided or generated separately.

- **Optional Flags:**  
  - Flags to enable or disable the calculation of partial volumes (when particles cross cell boundaries).  
  - A verbosity flag (`verbose`) can be set to provide diagnostic output during computation.

**Returns:**  
A tuple containing:
- `volume_per_cell_1` and `volume_per_cell_2`: Arrays of computed particle volume contributions for each cell in each dataset.

**Example:**
```
using Packing3D

volumes_per_cell_1, volumes_per_cell_2 = compute_volumes_per_cell(data_1, data_2, mesh=mesh)
```

---

### `get_cell_index`

**Description:**  
Determines the global cell index for a given coordinate vector using the provided mesh. The coordinate vector should match the expected format for the specified system (Cartesian: `[x, y, z]`; Cylindrical: `[r, theta, z]`).

**Positional and Keyword Arguments:**

- **Coordinate Input:**  
  - `coordinates::AbstractVector{Float64}` (Positional): A 3-element vector representing the position.

- **Mesh Input:**  
  - `mesh::Mesh` (keyword): The mesh object that defines the cell boundaries and divisions.

**Returns:**  
An integer representing the global index of the cell in which the coordinate is located.

**Example:**
```
using Packing3D

cell_idx = get_cell_index([0.5, 0.5, 0.5], mesh=mesh)
println("Cell Index: ", cell_idx)
```

---

### Mesh Utilities

The following types and functions help create and query the computational mesh:

- **`Mesh`**: A struct (with a constructor) representing the mesh (Cartesian or cylindrical).  
- **`get_mesh_boundaries`**: Returns an array of all cell boundaries in the mesh.  
- **`get_total_cells`**: Returns the total number of cells in the mesh.  
- **`get_cell_volume`**: Returns the volume of an individual cell.

**Example:**

```
using Packing3D

divisions = Dict(
  :r => 5,
  :theta => 3,
  :z => 10
)

params = Dict(
  :cylinder_radius => 0.03,
  :cylinder_base_level => 0.0,
  :cylinder_height => 0.08
)

mesh = Mesh(:cylindrical, divisions; params=params)

println("Total mesh cells: ", get_total_cells(mesh))
println("Cell volume in mesh: ", get_cell_volume(mesh))
```

---

### Data Partitioning Tools

- **`split_data`**, **`match_split_data`**, and **`extract_points`**:  
  Allow partitioning of a dataset into subsets based on spatial criteria or particle properties.

---

### Coordinate Conversion & I/O

- **`retrieve_coordinates`**: Extracts x, y, z coordinates and radii from a VTK data dictionary.  
- **`convert_to_cylindrical`**: Converts Cartesian coordinates to cylindrical coordinates for either scalars or arrays. 
  
  **Example:**

  ```
  using Packing3D

  data = read_vtk_file("particles.vtk")

  x_data, y_data, z_data, radii = retrieve_coordinates(data)

  r_data, theta_data = convert_to_cylindrical(x_data, y_data)
  ```

- **`read_vtk_file` and `save_vtk_file`**:  
  Provide robust reading from and writing to legacy ASCII POLYDATA VTK files.

  **Example:**

  ```
  using Packing3D

  data = read_vtk_file("particles.vtk")

  x_data, y_data, z_data, radii = retrieve_coordinates(data)

  params = Dict(
    :cylinder_radius => 0.03,
    :cylinder_base_level => 0.0,
    :cylinder_height => 0.08
  )

  mesh = Mesh(:cylindrical; params=params, divisions=Dict(:x => 10, :y => 10, :z => 10))
  
  total_volume_per_cell = compute_total_volume_per_cell(data; mesh=mesh)

  cell_volume = get_cell_volume(mesh)

  packing_densities = total_volume_per_cell ./ cell_volume

  num_particles = length(x_data)

  packing_densities = zeros(Float64, num_particles)

  for i in 1:num_particles
    cell_idx = get_cell_index([r_data[i], theta_data[i], z_data[i]]; mesh=mesh)

    packing_densities[i] = packing_density_per_cell[cell_idx]
  end
  
  data[:point_data][:packing_density] = packing_densities

  file_name = "particles_with_packing.vtk"
  save_vtk_file(file_name, data)
  println("Saved file: $file_name")

  ```

---

## How It Works

1. **Data Input:**  
   The package uses `read_vtk_file` to import particle data from legacy VTK files, extracting spatial coordinates and associated data.

2. **Coordinate Conversion:**  
   Use `convert_to_cylindrical` to transform Cartesian coordinates to cylindrical if needed.

3. **Meshing:**  
   The `Mesh` struct divides the defined region into cells for localised density calculations.

4. **Overlap Calculation:**  
   Analytical methods compute overlapping volumes for particles that span cell boundaries.

5. **Density & Mixing Metrics:**  
   Packing density is computed as the ratio of particle volume (including partial contributions) to the cell volume. Segregation and Lacey indices quantify the distribution of particle types.

---

## Requirements

- **Julia 1.6 or later**

Packing3D.jl is written entirely in Julia and uses only standard libraries from the Julia ecosystem.

---

## Installation

The package is not yet registered with Julia's package manager. To install, run:

```
using Pkg
Pkg.develop(url="https://github.com/fjbarter/packing3d_jl")
```

---

## Examples

Sample scripts demonstrating usage are provided in the `examples` directory. These include:
- Calculation of packing density in Cartesian and cylindrical systems.
- Splitting datasets based on spatial or particle property criteria.
- Computing segregation and Lacey mixing indices.

---

## Limitations

- **Particle Shape:**  
  Optimized for spherical particles.

- **VTK Format:**  
  Only legacy ASCII POLYDATA VTK files are supported. For XML VTK files, consider using [ReadVTK.jl](https://github.com/JuliaVTK/ReadVTK.jl).

---

## Planned Features

- **Enhanced Visualisation:**  
  Integration with interactive 3D plotting libraries.

- **Extended Particle Types:**  
  Future support for multisphere and non-spherical particles.

---

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

---

## Contact

For questions, support, or contributions, please contact:

**Freddie Barter**  
Email: [fjbarter@outlook.com](mailto:fjbarter@outlook.com)
