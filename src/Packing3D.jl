# Packing3D.jl # CASE SENSITIVE

# This is the parent module for the Packing3D package
module Packing3D

# Include mesh first
include("mesh.jl")
using .MeshModule
export Mesh, get_mesh_boundaries, get_total_cells, get_cell_volume

# Include each module
include("geometry.jl")
include("io.jl")
include("utils.jl") 
include("cartesian.jl")
include("cylindrical.jl")

# Bring in submodules (local imports)
using .Geometry
using .IO
using .Utils
using .Cartesian
using .Cylindrical

# Exports for public use
export calculate_packing,
       calculate_segregation_intensity,
       calculate_lacey,
       compute_volumes_per_cell,
       compute_total_volume_per_cell,
       get_cell_index,
       get_mesh_bounds,
       convert_to_cylindrical,
       read_vtk_file,
       save_vtk_file,
       split_data,
       match_split_data,
       retrieve_coordinates,
       extract_points
       

# Master functions for Lacey Index, Segregation Intensity, Packing Density

"""
    calculate_lacey(data_1, data_2; mesh::Union{Mesh, Nothing}=nothing, system::Symbol=:cartesian, 
                    params::Dict=Dict(), target_num_cells::Union{<:Real, Nothing}=1000, 
                    output_num_cells::Union{Bool, Nothing}=false, 
                    calculate_partial_volumes::Union{Bool, Nothing}=true, 
                    clamp_0_to_1::Union{Bool, Nothing}=true, verbose::Union{Bool, Nothing}=false) :: Float64

Compute the Lacey mixing index for two particle datasets using a mesh defined in either Cartesian or Cylindrical coordinates.
The function dispatches to the appropriate internal routine based on the chosen coordinate system.

# Arguments
- `data_1::Dict`: The first particle dataset.
- `data_2::Dict`: The second particle dataset.
- `mesh::Union{Mesh, Nothing}`: An optional preconstructed mesh. If not provided, one is generated using `target_num_cells` and `params`.
- `system::Symbol`: Coordinate system to use; either `:cartesian` or `:cylindrical` (default is `:cartesian`).
- `params::Dict`: Dictionary of mesh parameters.  
    - For Cartesian: must include `:x_min`, `:x_max`, `:y_min`, `:y_max`, `:z_min`, and `:z_max`.  
    - For Cylindrical: must include `:cylinder_radius`, `:cylinder_base_level`, and `:cylinder_height`.
- `target_num_cells`: Approximate desired number of mesh cells (default 1000).
- `output_num_cells`: If true, prints details on the target versus actual cell count.
- `calculate_partial_volumes`: If true, computes partial volumes for particles overlapping cell boundaries.
- `clamp_0_to_1`: If true, clamps the resulting index to the [0,1] range.
- `verbose`: If true, prints additional diagnostic output.

# Returns
- A `Float64` value representing the Lacey mixing index, where 0 indicates perfect mixing and 1 indicates complete segregation.

# Raises
- `ArgumentError` if required mesh parameters are missing.
"""
function calculate_lacey(
        data_1, data_2
        ;
        mesh::Union{Mesh, Nothing}=nothing,
        system::Union{Symbol, Nothing}=nothing,
        params::Dict=Dict(),
        target_num_cells::Union{<:Real, Nothing}=1000,
        output_num_cells::Union{Bool, Nothing}=false,
        calculate_partial_volumes::Union{Bool, Nothing}=true,
        clamp_0_to_1::Union{Bool, Nothing}=true,
        verbose::Union{Bool, Nothing}=false
    )::Float64

    if isnothing(mesh)
        if isnothing(target_num_cells)
            throw(ArgumentError("Target number of cells must be provided"))
        elseif isnothing(system)
            throw(ArgumentError("system must be provided if no mesh given. Choose :cartesian or :cylindrical"))
        else
            mesh = Mesh(system; target_num_cells=target_num_cells, params=params)
        end
    else
        system = mesh.system
    end

    if output_num_cells println("Target # of cells: $target_num_cells, Actual #: $(get_total_cells(mesh)), Divisions: $(mesh.divisions)") end

    if system == :cartesian
        
        lacey_index = _calculate_lacey_cartesian(
            data_1,
            data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            clamp_0_to_1=clamp_0_to_1,
            verbose=verbose
        )

    elseif system == :cylindrical

        lacey_index = _calculate_lacey_cylindrical(
            data_1,
            data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            clamp_0_to_1=clamp_0_to_1,
            verbose=verbose
        )

    end

    return lacey_index
end


"""
    calculate_segregation_intensity(data_1, data_2; mesh::Union{Mesh, Nothing}=nothing, system::Symbol=:cartesian, 
                                      params::Dict=Dict(), target_num_cells::Union{<:Real, Nothing}=1000, 
                                      output_num_cells::Union{Bool, Nothing}=false, 
                                      calculate_partial_volumes::Union{Bool, Nothing}=true, 
                                      clamp_0_to_1::Union{Bool, Nothing}=true, verbose::Union{Bool, Nothing}=false) :: Float64

Calculate the segregation intensity for two particle datasets using a mesh defined in Cartesian or Cylindrical coordinates.
This function delegates to the appropriate internal routine based on the selected coordinate system.

# Arguments
- `data_1::Dict`: Particle dataset for the first group.
- `data_2::Dict`: Particle dataset for the second group.
- `mesh::Union{Mesh, Nothing}`: An optional mesh structure. If not provided, a mesh is constructed using `target_num_cells` and `params`.
- `system::Symbol`: The coordinate system to use; must be either `:cartesian` or `:cylindrical` (default: `:cartesian`).
- `params::Dict`: Mesh parameters required for the chosen system.
- `target_num_cells`: Target number of mesh cells (default is 1000).
- `output_num_cells`: If true, prints mesh cell details.
- `calculate_partial_volumes`: If true, computes partial volumes for particles that span multiple cells.
- `clamp_0_to_1`: If true, clamps the resulting intensity to the range [0,1].
- `verbose`: If true, enables verbose diagnostic output.

# Returns
- A `Float64` value representing the segregation intensity (0 to 1).

# Raises
- `ArgumentError` if required parameters are missing.
"""
function calculate_segregation_intensity(
        data_1, data_2
        ;
        mesh::Union{Mesh, Nothing}=nothing,
        system::Symbol=:cartesian,
        params::Dict=Dict(),
        target_num_cells::Union{<:Real, Nothing}=1000,
        output_num_cells::Union{Bool, Nothing}=false,
        calculate_partial_volumes::Union{Bool, Nothing}=true,
        clamp_0_to_1::Union{Bool, Nothing}=true,
        verbose::Union{Bool, Nothing}=false
    )::Float64

    if isnothing(mesh)
        if isnothing(target_num_cells)
            throw(ArgumentError("Target number of cells must be provided"))
        else
            mesh = Mesh(system; target_num_cells=target_num_cells, params=params)
        end
    end

    if output_num_cells println("Target # of cells: $target_num_cells, Actual #: $(get_total_cells(mesh)), Divisions: $(mesh.divisions)") end

    if system == :cartesian
        
        segregation_intensity = _calculate_segregation_intensity_cartesian(
            data_1,
            data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            clamp_0_to_1=clamp_0_to_1,
            verbose=verbose
        )

    elseif system == :cylindrical

        segregation_intensity = _calculate_segregation_intensity_cylindrical(
            data_1,
            data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            clamp_0_to_1=clamp_0_to_1,
            verbose=verbose
        )

    end

    return segregation_intensity
end


"""
    calculate_packing(; file::Union{String, Nothing}=nothing, data::Union{Dict, Nothing}=nothing,
                      boundaries, system::Symbol = :cartesian,
                      cylinder_radius::Union{Float64, Nothing}=nothing,
                      accurate_cylindrical::Bool = true) -> Float64

Compute the packing density of particles within a defined region using either Cartesian or Cylindrical coordinates.
This function accepts either a particle data dictionary (via `data`) or a file path (via `file`) from which the particle
data is loaded. It also requires region boundaries and a coordinate system. If boundaries are provided as a dictionary, the
function verifies that the required keys are present and converts them to a vector using `convert_boundaries_dictionary`.

# Keyword Arguments
- `file::Union{String, Nothing}`: A path to a VTK file containing particle data. Used if `data` is not provided.
- `data::Union{Dict, Nothing}`: A dictionary containing particle data.
    - For Cartesian systems, expected keys are `:x_data`, `:y_data`, `:z_data`, and `:radii`.
    - For Cylindrical systems, expected keys are `:r_data`, `:theta_data`, `:z_data`, and `:radii`.
- `boundaries`: Either a `Dict{Symbol, <:Real}` or a `Vector{Float64}` defining the region boundaries.
    - For Cartesian systems, required keys (if a dictionary) are `:x_min`, `:x_max`, `:y_min`, `:y_max`, `:z_min`, and `:z_max`.
    - For Cylindrical systems, required keys (if a dictionary) are `:r_min`, `:r_max`, `:theta_min`, `:theta_max`, `:z_min`, and `:z_max`.
- `system::Symbol`: The coordinate system to use; must be either `:cartesian` (default) or `:cylindrical`.
- `cylinder_radius::Union{Float64, Nothing}`: (Optional) For Cartesian packing, if a cylindrical region is being considered.
- `accurate_cylindrical::Bool`: If true, uses an accurate cylindrical overlap calculation (default is true).

# Returns
- A `Float64` representing the packing density, defined as the fraction of the region's volume occupied by particles.

# Raises
- `ArgumentError` if neither `data` nor `file` is provided, if the provided boundaries (when given as a dictionary) are missing any required keys,
  or if the system is not recognized.
"""
function calculate_packing(; file::Union{String, Nothing}=nothing, data::Union{Dict, Nothing}=nothing,
                           boundaries, system::Symbol = :cartesian,
                           cylinder_radius::Union{Float64, Nothing}=nothing,
                           accurate_cylindrical::Bool = true,
                           centre::Tuple{<:Real, <:Real} = (0.0, 0.0),
                           calculate_partial_volumes::Bool = true) :: Float64
    # If no data is provided, attempt to load it from file.
    if isnothing(data)
        if isnothing(file)
            throw(ArgumentError("Either a file path or particle data must be provided."))
        else
            data = read_vtk_file(file)
        end
    end

    # Validate boundaries: if boundaries is a dictionary, check required keys.
    if isa(boundaries, Dict)
        if system == :cartesian
            required_keys = [:x_min, :x_max, :y_min, :y_max, :z_min, :z_max]
        elseif system == :cylindrical
            required_keys = [:r_min, :r_max, :theta_min, :theta_max, :z_min, :z_max]
        else
            throw(ArgumentError("Invalid system: $system. Must be :cartesian or :cylindrical."))
        end
        for key in required_keys
            if !haskey(boundaries, key)
                throw(ArgumentError("Boundaries missing required key: $key for $system system."))
            end
        end
        boundaries_vec = convert_boundaries_dictionary(boundaries, system)
    elseif isa(boundaries, Vector{Float64})
        boundaries_vec = boundaries
    else
        throw(ArgumentError("Invalid boundaries type; must be a Dict or Vector{Float64}."))
    end

    x_data, y_data, z_data, radii = retrieve_coordinates(data)

    if system == :cartesian
        # For Cartesian systems, data must contain :x_data, :y_data, :z_data, and :radii.
        packing_density = _compute_packing_cartesian(
            boundaries = boundaries_vec,
            x_data = x_data,
            y_data = y_data,
            z_data = z_data,
            radii = radii,
            cylinder_radius = cylinder_radius,
            calculate_partial_volumes = calculate_partial_volumes
        )
    elseif system == :cylindrical
        # For Cylindrical systems, data must contain :r_data, :theta_data, :z_data, and :radii.
        r_data, theta_data = convert_to_cylindrical(x_data, y_data; centre=centre)
        packing_density = _compute_packing_cylindrical(
            boundaries = boundaries_vec,
            r_data = r_data,
            theta_data = theta_data,
            z_data = z_data,
            radii = radii,
            accurate_cylindrical = accurate_cylindrical,
            calculate_partial_volumes = calculate_partial_volumes
        )
    else
        throw(ArgumentError("Invalid system: $system."))
    end

    return packing_density
end


"""
    compute_volumes_per_cell(data_1, data_2; mesh::Union{Mesh, Nothing}=nothing, 
                               calculate_partial_volumes::Union{Bool, Nothing}=true, verbose::Bool=false)

Compute the volume contributions for two particle datasets on a per-cell basis using the provided mesh.
This function returns the computed volume per cell for each dataset, along with the full particle volumes.

# Arguments
- `data_1::Dict`: The first particle dataset.
- `data_2::Dict`: The second particle dataset.
- `mesh::Union{Mesh, Nothing}`: A mesh structure specifying cell boundaries and divisions. Must be provided.
- `calculate_partial_volumes`: If true, calculates partial volumes for particles overlapping cell boundaries (default: true).
- `verbose`: If true, prints diagnostic output (default: false).

# Returns
- A tuple `(volume_per_cell_1, volume_per_cell_2, particle_volumes_1, particle_volumes_2)` where each element is an array of `Float64`.

# Raises
- `ArgumentError` if no mesh is provided.
"""
function compute_volumes_per_cell(
        data_1, data_2;
        mesh::Union{Mesh, Nothing}=nothing,
        calculate_partial_volumes::Union{Bool, Nothing}=true,
        verbose::Bool=false
    )

    if isnothing(mesh)
        throw(ArgumentError("Mesh must be provided"))
    end

    system = mesh.system

    if system == :cartesian

        volume_per_cell_1, volume_per_cell_2, real_pv_1, real_pv_2 = _compute_volume_per_cell_cartesian(
            data_1,
            data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            verbose=verbose
        )

        return volume_per_cell_1, volume_per_cell_2
        
    elseif system == :cylindrical

        volume_per_cell_1, volume_per_cell_2, real_pv_1, real_pv_2 = _compute_volume_per_cell_cylindrical(
            data_1,
            data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            verbose=verbose
        )

        return volume_per_cell_1, volume_per_cell_2

    else
        allowed_systems = [:cartesian, :cylindrical]
        throw(ArgumentError("Invalid coordinate system specified. Choose from: $allowed_systems"))
    end

end


"""
    compute_total_volume_per_cell(data; mesh::Union{Mesh, Nothing}=nothing, 
                                  calculate_partial_volumes::Union{Bool, Nothing}=true, verbose::Bool=false)

Compute the volume contributions for two particle datasets on a per-cell basis using the provided mesh.
This function returns the computed volume per cell for each dataset, along with the full particle volumes.

# Arguments
- `data::Dict`: The particle dataset.
- `mesh::Union{Mesh, Nothing}`: A mesh structure specifying cell boundaries and divisions. Must be provided.
- `calculate_partial_volumes`: If true, calculates partial volumes for particles overlapping cell boundaries (default: true).
- `verbose`: If true, prints diagnostic output (default: false).

# Returns
- An array `total_volume_per_cell` where each element is an array of `Float64`.

# Raises
- `ArgumentError` if no mesh is provided.
"""
function compute_total_volume_per_cell(
        data;
        mesh::Union{Mesh, Nothing}=nothing,
        calculate_partial_volumes::Union{Bool, Nothing}=true,
        verbose::Bool=false
    )

    if isnothing(mesh)
        throw(ArgumentError("Mesh must be provided"))
    end

    system = mesh.system

    num_points = length(data[:points][:, 1])
    mask = falses(num_points)
    mask[1] = true
    empty_data_2 = extract_points(data, mask)

    if system == :cartesian

        total_volume_per_cell, _, _, _ = _compute_volume_per_cell_cartesian(
            data,
            empty_data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            verbose=verbose
        )

        return total_volume_per_cell
        
    elseif system == :cylindrical

        total_volume_per_cell, _, _, _ = _compute_volume_per_cell_cylindrical(
            data,
            empty_data_2;
            mesh=mesh,
            calculate_partial_volumes=calculate_partial_volumes,
            verbose=verbose
        )

        return total_volume_per_cell

    else
        allowed_systems = [:cartesian, :cylindrical]
        throw(ArgumentError("Invalid coordinate system specified. Choose from: $allowed_systems"))
    end

end


"""
    get_cell_index(coordinates::AbstractVector{Float64}; mesh::Union{Mesh, Nothing}=nothing) -> Int

Determine the global cell index for a given coordinate vector based on the provided mesh.
This function computes the cell index using Cartesian or Cylindrical logic depending on the mesh system.

# Arguments
- `coordinates::AbstractVector{Float64}`: A vector of three coordinates.  
    - For Cartesian systems, these are `[x, y, z]`.  
    - For Cylindrical systems, these are `[r, theta, z]`.
- `mesh::Union{Mesh, Nothing}`: The mesh defining cell boundaries and divisions. Must be provided.

# Returns
- An `Int` representing the global index of the cell in which the coordinate lies.

# Raises
- `ArgumentError` if no mesh is provided or if the system is not recognized.
"""
function get_cell_index(coordinates::AbstractVector{Float64};
                        mesh::Union{Mesh, Nothing} = nothing)
    if mesh === nothing
        throw(ArgumentError("A mesh must be provided."))
    end

    system = mesh.system

    if system == :cartesian
        # Expect coordinates to be [x, y, z]
        x, y, z = coordinates[1:3]
        # Unpack parameters from the mesh.
        p = (; mesh.params...)  # Now p.x_min, p.x_max, etc.
        # Get the divisions (assumed stored with keys :x, :y, :z)
        d = mesh.divisions
        x_div, y_div, z_div = d[:x], d[:y], d[:z]
        # Compute cell sizes.
        dx = (p.x_max - p.x_min) / x_div
        dy = (p.y_max - p.y_min) / y_div
        dz = (p.z_max - p.z_min) / z_div
        recip_dx, recip_dy, recip_dz = 1/dx, 1/dy, 1/dz

        # Use the helper function (assumed available in Cartesian module)
        x_idx, y_idx, z_idx = Cartesian._compute_cell_index_cartesian(x, y, z, p.x_min, p.y_min, p.z_min, recip_dx, recip_dy, recip_dz)
        # Compute global index using another helper.
        global_idx = Cartesian._find_global_cell_index_cartesian(x_idx, y_idx, z_idx, x_div, y_div)
        # return (x_idx, y_idx, z_idx, global_idx)
        return global_idx
    elseif system == :cylindrical
        # Expect coordinates to be [r, theta, z]
        r, theta, z = coordinates[1:3]
        # Unpack cylindrical parameters from the mesh.
        p = (; mesh.params...)
        r_div = mesh.divisions[:r]
        z_div = mesh.divisions[:z]
        # For cylindrical, assume that p contains :cylinder_radius, :cylinder_base_level, :cylinder_height.
        cylinder_radius = p.cylinder_radius
        cylinder_base_level = p.cylinder_base_level
        cylinder_height = p.cylinder_height
        # Compute cell sizes.
        dr = cylinder_radius / r_div
        dz = cylinder_height / z_div

        # Use the helper function from the Cylindrical module.
        r_idx, theta_idx, z_idx = Cylindrical.compute_cell_index_cylindrical(r, theta, z, dr, r_div, dz, z_div, cylinder_base_level)
        global_idx = Cylindrical.find_global_cell_index_cylindrical([r_idx, theta_idx, z_idx], r_div, z_div)
        # return (r_idx, theta_idx, z_idx, global_idx)
        return global_idx
    else
        throw(ArgumentError("Invalid system: $system. Must be :cartesian or :cylindrical."))
    end
end

end # module Packing3D