module MeshModule

export Mesh, get_mesh_boundaries, get_total_cells, get_cell_boundaries, get_cell_volume, compute_divisions

"""
    Mesh(system::Symbol, divisions::Dict{String, Int}; params::Dict{String, Any}=Dict())

A structure representing a mesh in either Cartesian or Cylindrical coordinate systems.

# Arguments
- `system::Symbol`: Coordinate system, either `:cartesian` or `:cylindrical`.
- `divisions::Dict{String, Int}`: Number of divisions for each relevant dimension.
- `params::Dict{String, Any}`: Additional parameters required for mesh generation.

# Fields
- `system::Symbol`: The coordinate system.
- `divisions::Dict{String, Int}`: Number of divisions per dimension.
- `params::Dict{String, Any}`: Additional parameters.
- `total_cells::Int`: Total number of cells in the mesh.
- `cell_indices::Array{Int, 2}`: Indices mapping for each cell.
- `cell_boundaries::Array{Float64, 2}`: Boundaries defining each cell.
"""
struct Mesh
    system::Symbol                     # Coordinate system (:cartesian or :cylindrical)
    divisions::Dict{Symbol, Int}       # Number of divisions for each dimension
    params::Dict                       # Additional parameters
    total_cells::Int                   # Total number of cells
    cell_indices::Array{Int, 2}        # Indices for each cell
    cell_boundaries::Array{Float64, 2} # Boundaries for each cell
    cell_volume::Float64               # Volume of each cell

    # Inner constructor
    function Mesh(system::Symbol; target_num_cells::Union{<:Real, Nothing}=nothing, params::Dict=Dict(), divisions::Dict{Symbol, Int}=Dict{Symbol, Int}())

        # Validate params input
        validate_params!(system, params)

        if isempty(divisions)
            if isnothing(target_num_cells)
                throw(ArgumentError("Either target number of cells must be provided or divisions defined explicitly"))
            else
                divisions = compute_divisions(system; params=params, target_num_cells=target_num_cells)
            end
        end

        # Validate divisions input
        validate_divisions!(system, divisions)
        
        # Compute total number of cells based on the system and divisions
        total_cells = compute_total_cells(system, divisions, params)
        
        # Preallocate arrays for cell indices and boundaries
        cell_indices = Array{Int, 2}(undef, total_cells, 3)
        cell_boundaries = Array{Float64, 2}(undef, total_cells, 6)
        
        # Generate the mesh cells
        cell_volume = generate!(system, divisions, params, cell_indices, cell_boundaries)
        
        # Create and return the immutable Mesh instance
        new(system, divisions, params, total_cells, cell_indices, cell_boundaries, cell_volume)
    end
end

# Function to determine divisions from provided target number of cells
function compute_divisions(
        system::Symbol
        ;
        params::Dict=Dict(),
        target_num_cells::Real=1000
    )
    if system == :cartesian
        # Determine mesh divisions
        p = (; params...)  # Convert dict to named tuple
        x_min, x_max, y_min, y_max, z_min, z_max = p.x_min, p.x_max, p.y_min, p.y_max, p.z_min, p.z_max
        L_x = x_max - x_min
        L_y = y_max - y_min
        L_z = z_max - z_min

        # Ideal cubic cell side length
        s = (L_x * L_y * L_z / target_num_cells)^(1/3)

        # Divide each axis to make cells as cubic as possible
        division_vals = max.(1, round.(Int, [L_x/s, L_y/s, L_z/s]))
        
        keys = [:x, :y, :z]
        divisions = Dict{Symbol, Int}(zip(keys, division_vals))

        return divisions
    elseif system == :cylindrical
        # Determine mesh divisions
        z_divisions = max(1, round(Int, target_num_cells^(1 / 3)))
        num_cells_slice = target_num_cells / z_divisions
        theta_divisions = 3
        r_divisions = max(1, round(Int, sqrt(num_cells_slice)))

        divisions = Dict{Symbol, Int}(
            :r => r_divisions,
            :theta => theta_divisions,
            :z => z_divisions
        )

        return divisions
    else
        throw(ArgumentError("Invalid system: $system. Choose from :cartesian or :cylindrical"))
    end

end

# Validate that the parameters dictionary contains all required keys.
function validate_params!(system::Symbol, params::Dict)
    if system == :cartesian
        required_params = [:x_min, :x_max, :y_min, :y_max, :z_min, :z_max]
        for param in required_params
            if !haskey(params, param)
                throw(ArgumentError("Missing parameter: $param for Cartesian system."))
            end
        end
        validate_cartesian_params!(params)
    elseif system == :cylindrical
        required_params = [:cylinder_radius, :cylinder_base_level, :cylinder_height]
        for param in required_params
            if !haskey(params, param)
                throw(ArgumentError("Missing parameter: $param for Cylindrical system."))
            end
        end
        validate_cylindrical_params!(params)
    else
        throw(ArgumentError("Invalid system Symbol. Choose :cartesian or :cylindrical."))
    end
end

# Validate that the divisions dictionary contains all required keys.
function validate_divisions!(system::Symbol, divisions::Dict{Symbol, Int})
    if system == :cartesian
        for key in [:x, :y, :z]
            if !haskey(divisions, key)
                throw(ArgumentError("Missing division key: $key for Cartesian system. Found divisions: $divisions"))
            end
        end
    elseif system == :cylindrical
        for key in [:r, :theta, :z]
            if !haskey(divisions, key)
                throw(ArgumentError("Missing division key: $key for Cylindrical system."))
            end
        end
        # FOR NOW: enforce theta_div = 3 (required for constant-volume cells in the current mesh type)
        if divisions[:theta] != 3
            throw(ArgumentError("For constant_volume=true, 'theta' divisions must equal 3."))
        end
    end
end

# Function to validate Cartesian-specific parameters
function validate_cartesian_params!(params::Dict)
    if params[:x_min] >= params[:x_max]
        throw(ArgumentError("'x_min' must be less than 'x_max'."))
    end
    if params[:y_min] >= params[:y_max]
        throw(ArgumentError("'y_min' must be less than 'y_max'."))
    end
    if params[:z_min] >= params[:z_max]
        throw(ArgumentError("'z_min' must be less than 'z_max'."))
    end
end

# Function to validate Cylindrical-specific parameters
function validate_cylindrical_params!(params::Dict)
    if params[:cylinder_radius] <= 0
        throw(ArgumentError("'cylinder_radius' must be a positive number."))
    end
    if params[:cylinder_height] <= 0
        throw(ArgumentError("'cylinder_height' must be a positive number."))
    end
end

# Function to compute total number of cells
function compute_total_cells(system::Symbol, divisions::Dict{Symbol, Int}, params::Dict)::Int
    if system == :cartesian
        return divisions[:x] * divisions[:y] * divisions[:z]
    elseif system == :cylindrical
        r_div = divisions[:r]
        theta_div = divisions[:theta]
        z_div = divisions[:z]
        constant_volume = get(params, :constant_volume, true)
        if constant_volume
            # Total cells: r_div^2 * z_div
            return r_div^2 * z_div
        else
            # Compute number of slice cells as per Python logic
            slice_cells = 1 + sum(Int(round(theta_div * (2n - 1) / 3)) for n in 2:r_div)
            return slice_cells * z_div
        end
    end
end

# Function to generate mesh cells
function generate!(system::Symbol, divisions::Dict{Symbol, Int}, params::Dict,
                  cell_indices::Array{Int, 2}, cell_boundaries::Array{Float64, 2})
    if system == :cartesian
        return generate_cartesian_mesh!(divisions, params, cell_indices, cell_boundaries)
    elseif system == :cylindrical
        return generate_cylindrical_mesh!(divisions, params, cell_indices, cell_boundaries)
    end
end

# Function to generate Cartesian mesh cells
function generate_cartesian_mesh!(divisions::Dict{Symbol, Int}, params::Dict,
                                  cell_indices::Array{Int, 2}, cell_boundaries::Array{Float64, 2})
    x_div, y_div, z_div = divisions[:x], divisions[:y], divisions[:y]
    x_min, x_max = params[:x_min], params[:x_max]
    y_min, y_max = params[:y_min], params[:y_max]
    z_min, z_max = params[:z_min], params[:z_max]

    cell_volume = (x_max-x_min)*(y_max-y_min)*(z_max-z_min) / (x_div * y_div * z_div)

    # Create boundary ranges for each dimension
    x_bounds = LinRange(x_min, x_max, x_div + 1)
    y_bounds = LinRange(y_min, y_max, y_div + 1)
    z_bounds = LinRange(z_min, z_max, z_div + 1)

    index = 1
    for i in 1:x_div, j in 1:y_div, k in 1:z_div
        cell_indices[index, :] = [i, j, k]
        cell_boundaries[index, :] = [
            x_bounds[i], x_bounds[i + 1],
            y_bounds[j], y_bounds[j + 1],
            z_bounds[k], z_bounds[k + 1]
        ]
        index += 1
    end

    return cell_volume
end

# Function to generate Cylindrical mesh cells, including inner cells
function generate_cylindrical_mesh!(divisions::Dict{Symbol, Int}, params::Dict,
                                    cell_indices::Array{Int, 2}, cell_boundaries::Array{Float64, 2})
    r_div, theta_div, z_div = divisions[:r], divisions[:theta], divisions[:z]
    radius = params[:cylinder_radius]
    base_level = params[:cylinder_base_level]
    height = params[:cylinder_height]
    constant_volume = get(params, :constant_volume, true)

    # Calculate radius_inner
    radius_inner = radius / r_div

    cell_volume = pi * radius_inner^2 * (height / z_div)

    # Define radial bounds from 0 to radius, inclusive, with r_div radial divisions
    radial_bounds = LinRange(0.0, radius, r_div + 1)

    # Define vertical bounds
    z_bounds = LinRange(base_level, base_level + height, z_div + 1)

    index = 1

    for k in 1:z_div
        z_min, z_max = z_bounds[k], z_bounds[k + 1]

        if constant_volume
            # Add the special inner cell for this z-layer
            cell_indices[index, :] = [0, 0, k]
            cell_boundaries[index, :] = [
                -radius_inner, radius_inner,  # Radial boundaries (covers the center)
                0.0, 2 * π,                    # Angular boundaries (full circle)
                z_min, z_max                   # Vertical boundaries
            ]
            index += 1

            # Add all regular cells for this z-layer
            for i in 2:r_div
                r_min = radial_bounds[i]
                r_max = radial_bounds[i + 1]
                current_theta_div = 2 * i - 1  # 1, 3, 5, ..., 2*r_div -1

                # Create angular boundaries
                theta_bounds = LinRange(0.0, 2 * π, current_theta_div + 1)

                for j in 1:current_theta_div
                    theta_min, theta_max = theta_bounds[j], theta_bounds[j + 1]
                    cell_indices[index, :] = [i, j, k]
                    cell_boundaries[index, :] = [
                        r_min, r_max,
                        theta_min, theta_max,
                        z_min, z_max
                    ]
                    index += 1
                end
            end
        else
            # Variable volume: varying theta divisions per radial layer
            for i in 1:r_div
                if i == 1
                    # Add the special inner cell
                    cell_indices[index, :] = [0, 0, k]
                    cell_boundaries[index, :] = [
                        -radius_inner, radius_inner,
                        0.0, 2 * π,
                        z_min, z_max
                    ]
                    index += 1
                else
                    r_min = radial_bounds[i]
                    r_max = radial_bounds[i + 1]
                    # Calculate number of theta divisions for this radial layer
                    current_theta_div = max(1, round(Int, theta_div * (2 * i - 1) / 3))
                    
                    # Create angular boundaries
                    theta_bounds = LinRange(0.0, 2 * π, current_theta_div + 1)

                    for j in 1:current_theta_div
                        theta_min, theta_max = theta_bounds[j], theta_bounds[j + 1]
                        cell_indices[index, :] = [i, j, k]
                        cell_boundaries[index, :] = [
                            r_min, r_max,
                            theta_min, theta_max,
                            z_min, z_max
                        ]
                        index += 1
                    end
                end
            end
        end
    end

    # Trim the preallocated arrays if necessary
    if index - 1 != size(cell_boundaries, 1)
        throw(ArgumentError("Mismatch between computed cells ($index - 1) and allocated cells ($(size(cell_boundaries, 1))). Ensure that 'compute_total_cells' accurately reflects the mesh generation logic."))
    end

    return cell_volume
end

"""
    get_mesh_boundaries(mesh::Mesh) -> Array{Float64, 2}

Retrieve all cell boundaries in the mesh.

# Arguments
- `mesh::Mesh`: The mesh object.

# Returns
- `Array{Float64, 2}`: A 2D array where each row contains the boundaries of a cell in the format:
  `[r_min, r_max, theta_min, theta_max, z_min, z_max]` for cylindrical,
  or `[x_min, x_max, y_min, y_max, z_min, z_max]` for cartesian.
"""
function get_mesh_boundaries(mesh::Mesh)
    return mesh.cell_boundaries
end

"""
    get_cell_boundaries(mesh::Mesh, index::Int) -> Vector{Float64}

Retrieve the boundaries for a specific cell.

# Arguments
- `mesh::Mesh`: The mesh object.
- `index::Int`: The index of the cell.

# Returns
- `Vector{Float64}`: A vector containing the boundaries of the specified cell.

# Throws
- `ArgumentError` if the index is out of bounds.
"""
function get_cell_boundaries(mesh::Mesh, index::Int)
    if index < 1 || index > mesh.total_cells
        throw(ArgumentError("Cell index out of bounds. Valid range: 1 to $(mesh.total_cells)."))
    end
    return mesh.cell_boundaries[index, :]
end

"""
    get_total_cells(mesh::Mesh) -> Int

Retrieve the total number of cells in the mesh.

# Arguments
- `mesh::Mesh`: The mesh object.

# Returns
- `Int`: Total number of cells.
"""
function get_total_cells(mesh::Mesh)
    return mesh.total_cells
end


"""
    get_cell_volume(mesh::Mesh) -> Float64

Retrieve the volume of cells in the mesh.

# Arguments
- `mesh::Mesh`: The mesh object.

# Returns
- `Float64`: Volume of each cell.
"""
function get_cell_volume(mesh::Mesh)
    return mesh.cell_volume
end


end # module MeshModule
