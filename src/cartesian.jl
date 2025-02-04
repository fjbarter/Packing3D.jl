# cartesian.jl
# This module contains all the necessary functions that are unique to
# calculating packing density in cartesian coordinates

module Cartesian

# Import functions and structs from relevant modules
include("io.jl")
include("geometry.jl")
include("utils.jl")
include("mesh.jl")

# Using the modules (if functions are exported)
# using .IO, .Geometry, .Utils, .Cartesian, .MeshModule

# Import specific functions from modules
using .IO: read_vtk_file, retrieve_coordinates

using .Geometry: compute_cell_volume,
                 single_cap_intersection,
                 double_cap_intersection,
                 triple_cap_intersection

using .Utils: compute_automatic_boundaries,
              calculate_overlaps,
              calculate_active_overlap_values,
              is_inside_boundaries,
              is_outside_boundaries,
              convert_boundaries_dictionary

using .MeshModule: Mesh


export compute_packing_cartesian


"""
    compute_packing_cartesian(; file::Union{String, Nothing}=nothing,
                               boundaries::Union{Vector{Float64}, Nothing}=nothing,
                               x_data::Union{Vector{Float64}, Nothing}=nothing,
                               y_data::Union{Vector{Float64}, Nothing}=nothing,
                               z_data::Union{Vector{Float64}, Nothing}=nothing,
                               radii::Union{Vector{Float64}, Nothing}=nothing,
                               cylinder_radius::Union{Float64, Nothing}=nothing)

Compute the packing density of particles within a defined Cartesian region.

# Argumentsa
- `file::Union{String, Nothing}`: Path to the VTK file containing particle data.
  Required if `x_data`, `y_data`, `z_data`, and `radii` are not provided.
- `boundaries::Union{Vector{Float64}, Nothing}`: Array defining the Cartesian region boundaries.
- `x_data, y_data, z_data::Union{Vector{Float64}, Nothing}`: Particle coordinates.
- `radii::Union{Vector{Float64}, Nothing}`: Radii of the particles.
- `cylinder_radius::Union{Float64, Nothing}`: Radius of a cylindrical region (optional).

# Returns
- `Float64`: The packing density as the fraction of the Cartesian volume occupied by particles.

# Raises
- `ArgumentError`: If neither `file` nor particle data is provided.
"""
function compute_packing_cartesian(; file::Union{String, Nothing}=nothing,
                                    boundaries::Union{Vector{<:Real}, Dict{Symbol, <:Real}, Nothing}=nothing,
                                    x_data::Union{Vector{Float64}, Nothing}=nothing,
                                    y_data::Union{Vector{Float64}, Nothing}=nothing,
                                    z_data::Union{Vector{Float64}, Nothing}=nothing,
                                    radii::Union{Vector{Float64}, Nothing}=nothing,
                                    cylinder_radius::Union{Float64, Nothing}=nothing)::Float64
    # Step 1: Load data from file if not provided
    if x_data === nothing || y_data === nothing || z_data === nothing || radii === nothing
        if file === nothing
            throw(ArgumentError("Either 'file' or particle data must be provided."))
        end
        data = read_vtk_file(file)
        x_data, y_data, z_data, radii = retrieve_coordinates(data)
    end

    # Step 2: Determine boundaries
    if boundaries === nothing
        boundaries = compute_automatic_boundaries(; x_data=x_data, y_data=y_data, z_data=z_data, system=:cartesian)
    elseif isa(boundaries, Dict)
        boundaries = convert_boundaries_dictionary(boundaries, :cartesian)
    elseif !isa(boundaries, Vector{Float64})
        throw(ArgumentError("Boundaries must be a `Vector{Float64}`, `Dict{String, Float64}`, or `nothing`."))
    end

    # Step 3: Calculate overlaps
    overlaps = calculate_overlaps(; x_data=x_data, y_data=y_data, z_data=z_data,
                                    r_data=nothing, theta_data=nothing,
                                    radii=radii, boundaries=boundaries,
                                    factor=nothing,
                                    system=:cartesian)

    # Step 4: Calculate active overlap values
    total_particles = length(radii)
    active_overlap_values = calculate_active_overlap_values(total_particles; 
                                                            x_data=x_data, 
                                                            y_data=y_data, 
                                                            z_data=z_data,
                                                            r_data=nothing,
                                                            theta_data=nothing,
                                                            boundaries=boundaries, 
                                                            overlaps=overlaps, 
                                                            system=:cartesian)

    # Step 5: Pre-compute full particle volumes
    full_particle_volumes = (4 / 3) * π .* (radii .^ 3)

    # Step 6: Get masks for particles fully inside or outside boundaries
    inside_mask = is_inside_boundaries(; x_data=x_data, y_data=y_data, z_data=z_data,
                                         r_data=nothing, theta_data=nothing,
                                         boundaries=boundaries, radii=radii,
                                         factor=nothing,
                                         system=:cartesian)

    outside_mask = is_outside_boundaries(; x_data=x_data, y_data=y_data, z_data=z_data,
                                           r_data=nothing, theta_data=nothing,
                                           boundaries=boundaries, radii=radii,
                                           factor=nothing,
                                           system=:cartesian)

    # Step 7: Initialize total particle volume
    total_particle_volume = sum(full_particle_volumes[inside_mask])

    # Step 8: Compute volume for particles neither inside nor outside
    neither_mask = .~(inside_mask .| outside_mask)
    indices_neither = findall(neither_mask)
    
    # Step 9: Check if neither_mask is empty (all particles in or out)
    if !isempty(indices_neither)
        total_particle_volume += sum(calculate_particle_volume(i, radii, active_overlap_values) for i in indices_neither)
    end

    # Step 10: Compute cell volume
    cell_volume = compute_cell_volume(; boundaries=boundaries, system=:cartesian, cylinder_radius=cylinder_radius)

    # Step 11: Calculate packing density
    packing_density = total_particle_volume / cell_volume

    return packing_density
end


function calculate_particle_volume(i::Int, radii::Vector{Float64}, active_overlap_values::Matrix{Float64})::Float64
    # Initialize partial volume to zero
    partial_volume = 0.0

    # Extract valid overlap distances (non-NaN values)
    overlap_values = [val for val in active_overlap_values[i, :] if !isnan(val)]

    # Determine the number of overlaps for this particle
    number_of_overlaps = length(overlap_values)

    # Compute partial volume based on the number of overlaps
    if number_of_overlaps == 1
        # Single overlap: Use single-cap intersection geometry
        partial_volume = single_cap_intersection(radii[i], overlap_values[1])
    elseif number_of_overlaps == 2
        # Double overlap: Use double-cap intersection geometry
        partial_volume = double_cap_intersection(radii[i], overlap_values[1], overlap_values[2])
    elseif number_of_overlaps == 3
        # Triple overlap: Use triple-cap intersection geometry
        partial_volume = triple_cap_intersection(radii[i], overlap_values[1], overlap_values[2], overlap_values[3])
    end

    return partial_volume
end


end # modue Cartesian