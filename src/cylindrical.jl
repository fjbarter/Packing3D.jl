# cylindrical.jl
# This module contains all the necessary functions that are unique to
# calculating packing density in cylindrical coordinates

module Cylindrical

# Relative imports
using ..IO:         read_vtk_file, retrieve_coordinates

using ..Geometry:   convert_to_cylindrical,
                    calculate_angular_overlap_factor,
                    compute_cell_volume,
                    single_cap_intersection,
                    double_cap_intersection,
                    triple_cap_intersection,
                    sphere_cylinder_intersection,
                    sphere_cylinder_plane_intersection,
                    angular_difference

using ..Utils:      compute_automatic_boundaries,
                    convert_boundaries_dictionary

using ..Cartesian:  _calculate_partial_volume_cartesian

using ..MeshModule: Mesh,
                    get_mesh_boundaries,
                    get_total_cells,
                    get_cell_boundaries,
                    compute_divisions

# Exports to Packing3D top-level
export              _compute_packing_cylindrical,
                    _calculate_segregation_intensity_cylindrical,
                    _calculate_lacey_cylindrical,
                    _compute_volume_per_cell_cylindrical


"""
    compute_packing_cylindrical(; file::Union{String, Nothing}=nothing,
                                boundaries::Union{Vector{Float64}, Nothing}=nothing,
                                r_data::Union{Vector{Float64}, Nothing}=nothing,
                                theta_data::Union{Vector{Float64}, Nothing}=nothing,
                                z_data::Union{Vector{Float64}, Nothing}=nothing,
                                radii::Union{Vector{Float64}, Nothing}=nothing,
                                accurate_cylindrical::Bool=false)

Compute the packing density of particles within a cylindrical region.

# Arguments
- `file::Union{String, Nothing}`: Path to the VTK file containing particle data.
  Required if `r_data`, `theta_data`, `z_data`, and `radii` are not provided.
- `boundaries::Union{Vector{Float64}, Nothing}`: Array defining the cylindrical
  region boundaries: `[r_min, r_max, theta_min, theta_max, z_min, z_max]`.
- `r_data, theta_data, z_data::Union{Vector{Float64}, Nothing}`: Particle coordinates.
- `radii::Union{Vector{Float64}, Nothing}`: Radii of the particles.
- `accurate_cylindrical::Bool`: If true, computes accurate cylindrical overlaps.
  Otherwise, approximates overlaps as planar.

# Returns
- `Float64`: The packing density as the fraction of the cylindrical volume
  occupied by particles.

# Raises
- `ArgumentError`: If neither `file` nor particle data is provided.
"""
function _compute_packing_cylindrical(; file::Union{String, Nothing}=nothing,
                                      boundaries::Union{Vector{<:Real}, Dict{Symbol, <:Real}, Nothing}=nothing,
                                      r_data::Union{Vector{Float64}, Nothing}=nothing,
                                      theta_data::Union{Vector{Float64}, Nothing}=nothing,
                                      z_data::Union{Vector{Float64}, Nothing}=nothing,
                                      radii::Union{Vector{Float64}, Nothing}=nothing,
                                      accurate_cylindrical::Union{Bool, Nothing}=nothing,
                                      calculate_partial_volumes::Bool = true)::Float64

    # Step 1: Load data if necessary
    if r_data === nothing || theta_data === nothing || z_data === nothing || radii === nothing
        if file === nothing
            throw(ArgumentError("Either 'file' or particle data must be provided."))
        end
        try
            data = read_vtk_file(file)
            x_data, y_data, z_data, radii = retrieve_coordinates(data)
            r_data, theta_data = convert_to_cylindrical(x_data, y_data)
        catch e
            throw(e)
        end
    end

    # Step 2: Compute boundaries if not provided
    if boundaries === nothing
        boundaries = compute_automatic_boundaries(; r_data=r_data, theta_data=theta_data,
                                                   z_data=z_data, system=:cylindrical)
    else
        boundaries = convert_boundaries_dictionary(boundaries, :cylindrical)
    end

    # Find number of particles and check that all lengths are equal
    N = length(radii)
    if length(r_data) != N || length(theta_data) != N || length(z_data) != N
        throw(ArgumentError("r_data, theta_data, z_data and radii must all have the same length."))
    end
    # now it’s safe to @inbounds over 1:N

    # Compute the cell volume
    cell_volume = compute_cell_volume(; boundaries=boundaries, system=:cylindrical)

    # Short-circuit for centre-counting
    if !calculate_partial_volumes
        # Predefine accumulator for radius cubed (for volume calculation later)
        centre_inside_radii_cubed_sum = 0.0

        @inbounds for i in 1:N
            r, theta, z, radius = r_data[i], theta_data[i], z_data[i], radii[i]
            if _centre_inside_cylindrical(
                    r, theta, z, boundaries
                )

                centre_inside_radii_cubed_sum += radius^3
            end
        end

        return (4/3) * pi * centre_inside_radii_cubed_sum / cell_volume
    end

    r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

    if r_min == 0.0 && abs(theta_max - theta_min - 2π) < 1e-8
        throw(ArgumentError("r_min = 0 with a full circular theta range is invalid. Use r_min < 0 for a full cylinder."))
    end

    full_cylindrical_cell = r_min < 0.0

    # Enable accurate cylindrical calculation if the full cell condition is met
    if isnothing(accurate_cylindrical)
        accurate_cylindrical = full_cylindrical_cell ? true : false
    end

    # Then, define the partial volumes function depending on accurate_cylindrical
    partial_volume_fn = accurate_cylindrical ? 
        # this closure wraps the cylindrical version
        (radius, r, all_overlap_values, overlaps) -> _calculate_partial_volume_cylindrical(radius, r, all_overlap_values, overlaps) :
        # this one wraps the Cartesian (ignores the extra `r` argument)
        (radius, r, all_overlap_values, overlaps) -> _calculate_partial_volume_cartesian(radius, all_overlap_values, overlaps)

    completely_inside_radii_cubed_sum = 0.0
    total_partial_volume = 0.0

    @inbounds for i in 1:N
        r, theta, z, radius = r_data[i], theta_data[i], z_data[i], radii[i]
        if _completely_inside_cylindrical(
                r, theta, z, radius, boundaries
            )

            completely_inside_radii_cubed_sum += radius^3

        elseif _does_overlap_cylindrical(
                r, theta, z, radius, boundaries
            )
            
            all_overlap_values, overlaps = _particle_overlaps_cylindrical(
                boundaries, r, theta, z, radius
            )

            total_partial_volume += partial_volume_fn(radius, r, all_overlap_values, overlaps)

        # Else, particle completely outside => skip it
        end
    end

    total_particle_volume = (4/3) * pi * completely_inside_radii_cubed_sum + total_partial_volume

    return total_particle_volume/cell_volume

end


@inline function _centre_inside_cylindrical(
        r::Real,
        theta::Real,
        z::Real,
        boundaries::AbstractVector{<:Real}
    )

    # Extract boundaries
    r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

    # Full cylinder: no r_min constraint
    if r_min < 0 || isapprox(angular_difference(theta_min, theta_max), 2*pi; atol=1e-8)
        return (
            (r <= r_max) &&
            (z >= z_min) &&
            (z <= z_max)
        )
    end

    # Handle angular constraints
    theta_min = mod(theta_min, 2π)
    theta_max = mod(theta_max, 2π)

    standard_range = theta_min <= theta_max
    theta_inside = ifelse(
        standard_range,
        (theta >= theta_min) && (theta <= theta_max),
        (theta >= theta_min) || (theta <= theta_max)
    )

    return (
        (r >= r_min) &&
        (r <= r_max) &&
        theta_inside &&
        (z >= z_min) &&
        (z <= z_max)
    )
end


@inline function _completely_inside_cylindrical(
        r::Real,
        theta::Real,
        z::Real,
        radius::Real,
        boundaries::AbstractVector{<:Real}
    )
    # Extract boundaries
    r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

    # Full cylinder: no r_min constraint
    if r_min < 0
        return (
            (r <= r_max - radius) &&
            (z >= z_min + radius) &&
            (z <= z_max - radius)
        )
    end

    # Full ring: no angular range constraints
    if theta_min == 0 && theta_max == 2π
        return (
            (r >= r_min .+ radius) &&
            (r <= r_max .- radius) &&
            (z >= z_min .+ radius) &&
            (z <= z_max .- radius)
        )
    end

    # Handle angular constraints
    factor = _calculate_angular_overlap_factor(r, radius)

    theta_min = mod(theta_min + factor, 2π)
    theta_max = mod(theta_max - factor, 2π)

    # Theta condition for periodic ranges
    standard_range = (theta_min <= theta_max)
    theta_inside = ifelse(
        standard_range,
        (theta >= theta_min) && (theta <= theta_max),
        (theta >= theta_min) || (theta <= theta_max)
    )

    return (
        (r >= r_min + radius) &&
        (r <= r_max - radius) &&
        theta_inside &&
        (z >= z_min + radius) &&
        (z <= z_max - radius)
    )
end


@inline function _does_overlap_cylindrical(
        r::Real,
        theta::Real,
        z::Real,
        radius::Real,
        boundaries::AbstractVector{<:Real}
    )
    # Extract boundaries
    r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

    # Full cylinder: no r_min constraint
    if r_min < 0
        return (
            (r <= r_max + radius) &&
            (z >= z_min - radius) &&
            (z <= z_max + radius)
        )
    end

    # Full ring: no angular range constraints
    if theta_min == 0 && theta_max == 2π
        return (
            (r >= r_min - radius) &&
            (r <= r_max + radius) &&
            (z >= z_min - radius) &&
            (z <= z_max + radius)
        )
    end

    factor = _calculate_angular_overlap_factor(r, radius)

    theta_min = mod(theta_min - factor, 2π)
    theta_max = mod(theta_max + factor, 2π)

    # Theta condition for periodic ranges
    standard_range = (theta_min <= theta_max)
    theta_inside = ifelse(
        standard_range,
        (theta >= theta_min) && (theta <= theta_max),
        (theta >= theta_min) || (theta <= theta_max)
    )

    return (
        (r >= r_min - radius) &&
        (r <= r_max + radius) &&
        theta_inside &&
        (z >= z_min - radius) &&
        (z <= z_max + radius)
    )
end


@inline function _calculate_angular_overlap_factor(r, radius)

    # Ensure radial distances are safe (prevent division by zero)
    safe_r = max(r, radius)

    # Calculate angular factor using arcsin, ensuring the ratio is clipped to [-1, 1]
    factor = asin(clamp(radius / safe_r, -1.0, 1.0))

    return factor
end


@inline function _calculate_partial_volume_cylindrical(
        radius::Float64,
        r::Float64,
        all_overlap_values::Vector{Float64},
        overlaps::BitVector,
    )::Float64

    # radial overlap?
    overlaps_r = overlaps[1] || overlaps[2]
    if overlaps[1]
        r_overlap, min_boundary_r = all_overlap_values[1], true
    elseif overlaps[2]
        r_overlap, min_boundary_r = all_overlap_values[2], false
    end

    # angular overlap?
    overlaps_theta = overlaps[3] || overlaps[4]
    if overlaps[3]
        theta_overlap = all_overlap_values[3]
    elseif overlaps[4]
        theta_overlap = all_overlap_values[4]
    end

    # axial (z) overlap?
    overlaps_z = overlaps[5] || overlaps[6]
    if overlaps[5]
        z_overlap = all_overlap_values[5]
    elseif overlaps[6]
        z_overlap = all_overlap_values[6]
    end

    # Count the number of overlaps
    num_overlaps = (overlaps_r ? 1 : 0) +
                   (overlaps_theta ? 1 : 0) +
                   (overlaps_z ? 1 : 0)

    # Calculate partial volume based on the number of overlaps
    if num_overlaps == 1
        if overlaps_r
            partial_volume = sphere_cylinder_intersection(
                radius, r, r_overlap; min_boundary=min_boundary_r
            )
        elseif overlaps_theta
            partial_volume = single_cap_intersection(radius, theta_overlap)
        elseif overlaps_z
            partial_volume = single_cap_intersection(radius, z_overlap)
        end
    elseif num_overlaps == 2
        if overlaps_r && overlaps_theta
            partial_volume = double_cap_intersection(radius, r_overlap, theta_overlap)
        elseif overlaps_r && overlaps_z
            partial_volume = sphere_cylinder_plane_intersection(
                radius, r, r_overlap, z_overlap; min_boundary=min_boundary_r
            )
        elseif overlaps_theta && overlaps_z
            partial_volume = double_cap_intersection(radius, theta_overlap, z_overlap)
        end
    elseif num_overlaps == 3
        if overlaps_r && overlaps_theta && overlaps_z
            partial_volume = triple_cap_intersection(
                radius, r_overlap, theta_overlap, z_overlap
            )
        end
    else
        # No valid overlaps
        partial_volume = 0.0
    end

    return partial_volume
end


"""
    calculate_segregation_intensity(data_1::Dict, data_2::Dict, cylinder_radius::Float64,
                                    cylinder_base_level::Float64, cylinder_height::Float64,
                                    target_num_cells::Int; packing_threshold::Float64=0.05)

Calculate the segregation intensity for two particle datasets within a cylindrical mesh.

# Arguments
- `data_1::Dict`: Particle dataset for the first group (e.g., small particles).
- `data_2::Dict`: Particle dataset for the second group (e.g., large particles).
- `cylinder_radius::Float64`: Radius of the cylindrical region.
- `cylinder_base_level::Float64`: Base level of the cylindrical region.
- `cylinder_height::Float64`: Height of the cylindrical region.
- `target_num_cells::Int`: Approximate number of cells in the cylindrical mesh.
- `packing_threshold::Float64`: Minimum packing density to consider a cell valid.

# Returns
- `Float64`: Segregation intensity, a dimensionless value ranging from 0 (perfectly mixed) to 1 (completely segregated). Returns `NaN` if no valid cells are found.
"""
function _calculate_segregation_intensity_cylindrical(data_1::Dict,
                                         data_2::Dict
                                         ;
                                         mesh::Union{Mesh, Nothing}=nothing,
                                         calculate_partial_volumes::Bool=true,
                                         clamp_0_to_1::Bool=true,
                                         verbose::Bool=false)::Float64

    volume_per_cell_1, volume_per_cell_2, real_pv_1, real_pv_2 = _compute_volume_per_cell_cylindrical(
        data_1,
        data_2;
        mesh=mesh,
        calculate_partial_volumes=calculate_partial_volumes,
        verbose=verbose
    )
    # Calculate concentration of type 1 in each cell (type agnostic)
    total_particle_volume_per_cell = (volume_per_cell_1 .+ volume_per_cell_2)

    total_particle_volume = sum(total_particle_volume_per_cell)

    min_particle_volume = minimum(vcat(real_pv_1, real_pv_2))

    if verbose
        println("Total of assigned cell volumes: $(sum(total_particle_volume_per_cell))")
        println("Actual total particle volume:   $(sum(real_pv_1) + sum(real_pv_2))")
    end

    problem_cells = findall(total_particle_volume_per_cell .< -min_particle_volume)
    if !isempty(problem_cells)
        # println("WARNING: Total particle volume in cell calculated as negative at indices: ", problem_cells)
        # println("Volume calculated as: ", total_particle_volume_per_cell[problem_cells])
        return NaN
    end


    if total_particle_volume == 0
        println("Total particle volume calculated as zero, returning NaN")
        return NaN
    end

    conc_mean = sum(volume_per_cell_1) / total_particle_volume

    # Calculate concentration for occupied cells (0 if unoccupied)
    concs_1 = ifelse.(
        total_particle_volume_per_cell .== 0,
        0.0,
        volume_per_cell_1 ./ total_particle_volume_per_cell
    )

    # Step 8: Calculate segregation intensity
    weighted_numerator = sum(total_particle_volume_per_cell .* ((concs_1 .- conc_mean).^2)) / total_particle_volume
    I_S_max = sqrt(conc_mean * (1 - conc_mean))
    segregation_intensity = sqrt(weighted_numerator) / I_S_max


    return segregation_intensity
end


"""
    calculate_lacey(data_1::Dict, data_2::Dict, cylinder_radius::Float64, 
                    cylinder_base_level::Float64, cylinder_height::Float64, 
                    target_num_cells::Int; packing_threshold::Float64=0.05,
                    output_num_cells::Bool=false)

Calculate the Lacey mixing index for two particle datasets in a cylindrical mesh.

# Arguments
- `data_1::Dict`: Particle dataset for the first group (e.g., small particles).
- `data_2::Dict`: Particle dataset for the second group (e.g., large particles).
- `cylinder_radius::Float64`: Radius of the cylindrical region.
- `cylinder_base_level::Float64`: Base level of the cylindrical region.
- `cylinder_height::Float64`: Height of the cylindrical region.
- `target_num_cells::Int`: Approximate number of cells in the cylindrical mesh.
- `packing_threshold::Float64`: Minimum packing density to consider a cell valid.
- `output_num_cells::Bool`: If true, prints the actual number of cells in the mesh.

# Returns
- `Float64`: The Lacey mixing index, ranging from 0 (perfect mixing) to 1 (complete segregation).
"""
function _calculate_lacey_cylindrical(data_1::Dict,
                         data_2::Dict
                         ;
                         mesh::Union{Mesh, Nothing}=nothing,
                         calculate_partial_volumes::Union{Bool, Nothing}=true,
                         clamp_0_to_1::Union{Bool, Nothing}=true,
                         verbose::Union{Bool, Nothing}=false)::Float64

    volume_per_cell_1, volume_per_cell_2, real_pv_1, real_pv_2 = _compute_volume_per_cell_cylindrical(
        data_1,
        data_2;
        mesh=mesh,
        calculate_partial_volumes=calculate_partial_volumes,
        verbose=verbose
    )
    
    # Calculate concentration of type 1 in each cell (type agnostic)
    total_particle_volume_per_cell = (volume_per_cell_1 .+ volume_per_cell_2)

    total_particle_volume = sum(total_particle_volume_per_cell)

    min_particle_volume = minimum(vcat(real_pv_1, real_pv_2))

    if verbose
        println("Total of assigned cell volumes: $(sum(total_particle_volume_per_cell))")
        println("Actual total particle volume:   $(sum(real_pv_1) + sum(real_pv_2))")
    end

    problem_cells = findall(total_particle_volume_per_cell .< -min_particle_volume)
    if !isempty(problem_cells)
        # println("WARNING: Total particle volume in cell calculated as negative at indices: ", problem_cells)
        # println("Volume calculated as: ", total_particle_volume_per_cell[problem_cells])
        return NaN
    end


    if total_particle_volume == 0
        println("Total particle volume calculated as zero, returning NaN")
        return NaN
    end

    conc_mean = sum(volume_per_cell_1) / total_particle_volume

    # Calculate concentration for occupied cells (0 if unoccupied)
    concs_1 = ifelse.(
        total_particle_volume_per_cell .== 0,
        0.0,
        volume_per_cell_1 ./ total_particle_volume_per_cell
    )

    # Step 7: Compute effective particles per cell
    # eff_particles_per_cell = sum(total_particle_volume_per_cell)^2 / sum(total_particle_volume_per_cell.^2)

    # println("lenconcs: $(length(concs_1_valid)), num_valid_cells: $(num_valid_cells)")

    # Step 10: Compute Lacey index
    # Actual variance of concentrations (weighted by cell occupanacy)
    S_actual = sum(total_particle_volume_per_cell .* ((concs_1 .- conc_mean).^2)) / total_particle_volume

    # Random variance (binomial distribution assumption)
    # S_random = conc_mean * (1 - conc_mean) / eff_particles_per_cell
    S_random = 0
    
    # Maximum variance (complete segregation, binomial distribution assumption)
    S_maximum = conc_mean * (1 - conc_mean)

    if verbose
        println("Total of assigned cell volumes: $(sum(total_particle_volume_per_cell))")
        println("Actual total particle volume:   $(sum(volume_per_cell_1) + sum(volume_per_cell_2))")
        println("S_random: $S_random")
        println("S_maximum: $S_maximum")
        println("S_actual: $S_actual")
    end

    # Final calculation of Lacey index
    M = (S_maximum - S_actual) / (S_maximum - S_random)

    if clamp_0_to_1
        M = clamp(M, 0, 1)
    end

    return M
end


function _compute_volume_per_cell_cylindrical(data_1::Dict,
                                 data_2::Dict
                                 ;
                                 mesh::Union{Mesh, Nothing}=nothing,
                                 calculate_partial_volumes::Union{Bool, Nothing}=true,
                                 verbose::Union{Bool, Nothing}=false)

    # Extract mesh Information
    params = mesh.params
    divisions = mesh.divisions
    if haskey(params, :centre)
        centre = params[:centre]
    else
        centre = (0.0, 0.0)
    end

    # Extract particle data from both datasets
    x_data_1, y_data_1, z_data_1, radii_1 = retrieve_coordinates(data_1)
    x_data_2, y_data_2, z_data_2, radii_2 = retrieve_coordinates(data_2)

    # Convert Cartesian to cylindrical coordinates
    r_data_1, theta_data_1 = convert_to_cylindrical(x_data_1, y_data_1; centre=centre)
    r_data_2, theta_data_2 = convert_to_cylindrical(x_data_2, y_data_2; centre=centre)

    # Step 3: Compute total volumes of both datasets
    particle_volumes_1 = (4 / 3) * π * (radii_1 .^ 3)
    particle_volumes_2 = (4 / 3) * π * (radii_2 .^ 3)

    max_particle_diameter = 2 * maximum(vcat(radii_1, radii_2))

    cylinder_radius, cylinder_base_level, cylinder_height = params[:cylinder_radius], params[:cylinder_base_level], params[:cylinder_height]
    r_divisions, z_divisions = divisions[:r], divisions[:z]
    
    radius_inner = cylinder_radius / r_divisions
    cell_volume = π * radius_inner^2 * cylinder_height / z_divisions

    cell_too_small = min(cylinder_radius/r_divisions, cylinder_height/z_divisions) < max_particle_diameter

    if calculate_partial_volumes && cell_too_small
        if verbose println("WARNING: Mesh resolution is too fine; refusing to calculate partial volumes") end
        calculate_partial_volumes = false
    end

    mesh_boundaries = get_mesh_boundaries(mesh)
    num_cells = get_total_cells(mesh)

    num_particles_1 = length(r_data_1)
    num_particles_2 = length(r_data_2)

    # Step 5: Initialise zeros array for cell volumes
    total_particle_volume_per_cell = fill(0.0, num_cells)

    # Calculate division size for finding cell index from positions
    dr = cylinder_radius / r_divisions
    dz = cylinder_height / z_divisions

    # Step 6: Decide whether to calculate partial volumes or not
    if calculate_partial_volumes
        # Calculate partial volume contribution to overlapped cells
        # Refer to partial volume calculation function
        volume_per_cell_1 = compute_partial_volume_per_cell_cylindrical(mesh_boundaries, num_cells, num_particles_1, particle_volumes_1, r_data_1, theta_data_1, z_data_1, radii_1, dr, r_divisions, dz, z_divisions, cylinder_base_level)
        volume_per_cell_2 = compute_partial_volume_per_cell_cylindrical(mesh_boundaries, num_cells, num_particles_2, particle_volumes_2, r_data_2, theta_data_2, z_data_2, radii_2, dr, r_divisions, dz, z_divisions, cylinder_base_level)
    else
        # Not calculating partial volumes
        # Bin particles based on their centre positions

        volume_per_cell_1 = compute_centres_volume_per_cell_cylindrical(num_cells, num_particles_1, particle_volumes_1, r_data_1, theta_data_1, z_data_1, dr, r_divisions, dz, z_divisions, cylinder_base_level)
        volume_per_cell_2 = compute_centres_volume_per_cell_cylindrical(num_cells, num_particles_2, particle_volumes_2, r_data_2, theta_data_2, z_data_2, dr, r_divisions, dz, z_divisions, cylinder_base_level)

    end

    return volume_per_cell_1, volume_per_cell_2, particle_volumes_1, particle_volumes_2
end


function compute_partial_volume_per_cell_cylindrical(mesh_boundaries, num_cells, num_particles, particle_volumes, r_data, theta_data, z_data, radii, dr, r_divisions, dz, z_divisions, cylinder_base_level)
    # Initialise particle_volume_per_cell
    particle_volume_per_cell = zeros(Float64, num_cells)

    # case_counter = zeros(Int, 9)
    # total_case_counter = zeros(Int, 9)
    
    # Loop through all cells and calculate packing densities, then concentration of species 1
    for i in 1:num_particles
        r_idx, theta_idx, z_idx = compute_cell_index_cylindrical(r_data[i], theta_data[i], z_data[i], dr, r_divisions, dz, z_divisions, cylinder_base_level)
        if _is_index_outside_cylindrical(r_idx, z_idx, r_divisions, z_divisions)
            continue
        end
        cell_idx = find_global_cell_index_cylindrical([r_idx, theta_idx, z_idx], r_divisions, z_divisions)
        overlap_values, overlaps = _particle_overlaps_cylindrical(mesh_boundaries[cell_idx, :], r_data[i], theta_data[i], z_data[i], radii[i])
        num_overlaps = sum(overlaps)

        overlaps_r = any(overlaps[1:2])
        overlaps_theta = any(overlaps[3:4])
        overlaps_z = any(overlaps[5:6])

        if overlaps[2] && (r_idx == r_divisions)
            # println("Warning: r_max overlap calculated in outermost division")
            # println("Coordinates: $(r_data[i]), $(theta_data[i]), $(z_data[i])")
            # Pretend the whole particle is in the cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        end

        if overlaps[5] && (z_idx == 1)
            # println("Warning: z_min overlap calculated in bottom division")
            # println("Coordinates: $(r_data[i]), $(theta_data[i]), $(z_data[i]), radius $(radii[i])")
            # println("Overlap values: $overlap_values")
            # println("Overlaps mask: $overlaps")
            # println("Mesh boundaries: $(mesh_boundaries[cell_idx, :])")
            # Pretend the whole particle is in the cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        end

        if overlaps[6] && (z_idx == z_divisions)
            # println("Warning: z_max overlap calculated in top division")
            # println("Coordinates: $(r_data[i]), $(theta_data[i]), $(z_data[i])")
            # Pretend the whole particle is in the cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        end

        relative_r, relative_theta, relative_z = 0, 0, 0

        if overlaps_r
            if overlaps[1]
                relative_r = -1
            else
                relative_r = 1
            end
            adjacent_overlap, adjacent_theta_idx, adjacent_relative_theta = adjacent_radial_overlaps(r_data[i], theta_data[i], radii[i], r_idx, relative_r)
            min_boundary = Bool((1 - relative_r) / 2)
        end
        if overlaps_theta
            if overlaps[3]
                relative_theta = -1
            else
                relative_theta = 1
            end
        end
        if overlaps_z
            if overlaps[5]
                relative_z = -1
            else
                relative_z = 1
            end
        end

        # case = 0



        if num_overlaps == 0
            # Particle completely inside cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        elseif num_overlaps == 1
            overlap_value = overlap_values[overlaps][1]
            if overlaps_r
                # main_volume = sphere_cylinder_intersection(radii[i], r_data[i], overlap_value; min_boundary=min_boundary)
                main_volume = single_cap_intersection(radii[i], overlap_value)
                if isnothing(adjacent_overlap)
                    remainder_volumes = [particle_volumes[i] - main_volume]
                    remainder_indices = [[r_idx + relative_r, adjacent_theta_idx, z_idx]]
                else
                    # case = 1
                    remainder_volumes = zeros(Float64, 2)
                    remainder_volumes[1] = double_cap_intersection(radii[i], -overlap_value, adjacent_overlap)
                    remainder_volumes[2] = particle_volumes[i] - main_volume - remainder_volumes[1]
                    remainder_indices = [[r_idx + relative_r, adjacent_theta_idx, z_idx],
                                         [r_idx + relative_r, adjacent_theta_idx + adjacent_relative_theta, z_idx]]
                end
            else
                # Overlaps z or theta boundary, giving a single spherical cap
                main_volume = single_cap_intersection(radii[i], overlap_value)
                remainder_volumes = [particle_volumes[i] - main_volume]

                remainder_indices = [[r_idx + relative_r, theta_idx + relative_theta, z_idx + relative_z]]
            end
        elseif num_overlaps == 2
            # Make assumption that curvature near corners is minimal
            # Therefore, double spherical caps in all four corners
            overlap_vals = overlap_values[overlaps]

            main_volume = double_cap_intersection(radii[i], overlap_vals[1], overlap_vals[2])
            
            # Three remainders, need to find their volumes

            if overlaps_theta && overlaps_z
                # case = 2
                remainder_indices = fill(zeros(Int, 3), 3)
                remainder_volumes = zeros(Float64, 3)
                # Three remainders
                remainder_indices[1] = [r_idx, theta_idx + relative_theta, z_idx]
                remainder_indices[2] = [r_idx, theta_idx, z_idx + relative_z]
                remainder_indices[3] = [r_idx, theta_idx + relative_theta, z_idx + relative_z]
                
                remainder_volumes[1] = single_cap_intersection(radii[i], overlap_vals[2]) - main_volume
                remainder_volumes[2] = single_cap_intersection(radii[i], overlap_vals[1]) - main_volume
                remainder_volumes[3] = particle_volumes[i] - remainder_volumes[1] - remainder_volumes[2] - main_volume
            elseif overlaps_r && overlaps_theta
                if isnothing(adjacent_overlap)
                    # case = 3
                    remainder_indices = fill(zeros(Int, 2), 2)
                    remainder_volumes = zeros(Float64, 2)
                    # remainder_volumes[1] = sphere_cylinder_intersection(radii[i], r_data[i], -overlap_vals[1]; min_boundary=!min_boundary)
                    remainder_volumes[1] = single_cap_intersection(radii[i], -overlap_vals[1])
                    remainder_indices[1] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                    remainder_volumes[2] = particle_volumes[i] - main_volume - remainder_volumes[1]
                    remainder_indices[2] = [r_idx, theta_idx + relative_theta, z_idx]
                else
                    # case = 4
                    remainder_indices = fill(zeros(Int, 3), 3)
                    remainder_volumes = zeros(Float64, 3)
                    remainder_volumes[1] = single_cap_intersection(radii[i], overlap_vals[1]) - main_volume
                    remainder_indices[1] = [r_idx, theta_idx + relative_theta, z_idx]
                    remainder_volumes[2] = double_cap_intersection(radii[i], -overlap_vals[1], adjacent_overlap)
                    remainder_indices[2] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                    remainder_volumes[3] = particle_volumes[i] - main_volume - remainder_volumes[1] - remainder_volumes[2]
                    remainder_indices[3] = [r_idx + relative_r, adjacent_theta_idx + adjacent_relative_theta, z_idx]
                end
                
            elseif overlaps_r && overlaps_z
                if isnothing(adjacent_overlap)
                    # case = 5
                    remainder_indices = fill(zeros(Int, 3), 3)
                    remainder_volumes = zeros(Float64, 3)
                    # remainder_volumes[1] = sphere_cylinder_intersection(radii[i], r_data[i], overlap_vals[1]) - main_volume
                    remainder_volumes[1] = single_cap_intersection(radii[i], overlap_vals[1]) - main_volume
                    remainder_indices[1] = [r_idx, theta_idx, z_idx + relative_z]
                    remainder_volumes[2] = double_cap_intersection(radii[i], -overlap_vals[1], overlap_vals[2])
                    remainder_indices[2] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                    remainder_volumes[3] = particle_volumes[i] - main_volume - remainder_volumes[1] - remainder_volumes[2]
                    remainder_indices[3] = [r_idx + relative_r, adjacent_theta_idx, z_idx + relative_z]
                else
                    # case = 6
                    remainder_indices = fill(zeros(Int, 5), 5)
                    remainder_volumes = zeros(Float64, 5)
                    # remainder_volumes[1] = sphere_cylinder_intersection(radii[i], r_data[i], overlap_vals[1]) - main_volume
                    remainder_volumes[1] = single_cap_intersection(radii[i], overlap_vals[1]) - main_volume
                    remainder_indices[1] = [r_idx, theta_idx, z_idx + relative_z]
                    remainder_volumes[2] = triple_cap_intersection(radii[i], -overlap_vals[1], overlap_vals[2], adjacent_overlap)
                    remainder_indices[2] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                    remainder_volumes[3] = triple_cap_intersection(radii[i], -overlap_vals[1], overlap_vals[2], -adjacent_overlap)
                    remainder_indices[3] = [r_idx + relative_r, adjacent_theta_idx + adjacent_relative_theta, z_idx]
                    remainder_volumes[4] = triple_cap_intersection(radii[i], -overlap_vals[1], -overlap_vals[2], adjacent_overlap)
                    remainder_indices[4] = [r_idx + relative_r, adjacent_theta_idx, z_idx + relative_z]
                    remainder_volumes[5] = particle_volumes[i] - main_volume - remainder_volumes[1] - remainder_volumes[2] - remainder_volumes[3] - remainder_volumes[4]
                    remainder_indices[5] = [r_idx + relative_r, adjacent_theta_idx + adjacent_relative_theta, z_idx + relative_z]
                end
            else
                # println("Calculated 2 overlaps but could not match the case")
                # println("Overlap values: $overlap_values")
                # println("Overlaps mask: $overlaps")
                # println("Mesh boundaries: $(mesh_boundaries[cell_idx, :])")
                # println("Coordinates: $(r_data[i]), $(theta_data[i]), $(z_data[i])")
                # Pretend its inside
                main_volume = particle_volumes[i]
                particle_volume_per_cell[cell_idx] += main_volume
                continue
            end
        elseif num_overlaps == 3
            # Not going to bother working out triple caps (contribution is so small)
            # Make assumption that average position is on the corner -> 1 eighth in each direction
            # Total volume still conserved
            main_volume = particle_volumes[i] / 8
            if isnothing(adjacent_overlap)
                # case = 7
                remainder_indices = fill(zeros(Int, 5), 5)
                remainder_volumes = zeros(Float64, 5)
                remainder_volumes[1:3] .= main_volume
                remainder_indices[1] = [r_idx, theta_idx, z_idx + relative_z]
                remainder_indices[2] = [r_idx, theta_idx + relative_theta, z_idx + relative_z]
                remainder_indices[3] = [r_idx, theta_idx + relative_theta, z_idx]
                remainder_volumes[4:5] .= main_volume*2
                remainder_indices[4] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                remainder_indices[5] = [r_idx + relative_r, adjacent_theta_idx, z_idx + relative_z]
            else
                # case = 8
                remainder_indices = fill(zeros(Int, 7), 7)
                remainder_volumes = zeros(Float64, 7)
                remainder_volumes[1:7] .= main_volume
                remainder_indices[1] = [r_idx, theta_idx, z_idx + relative_z]
                remainder_indices[2] = [r_idx, theta_idx + relative_theta, z_idx]
                remainder_indices[3] = [r_idx, theta_idx + relative_theta, z_idx + relative_z]
                remainder_indices[4] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                remainder_indices[5] = [r_idx + relative_r, adjacent_theta_idx + adjacent_relative_theta, z_idx]
                remainder_indices[6] = [r_idx + relative_r, adjacent_theta_idx, z_idx + relative_z]
                remainder_indices[7] = [r_idx + relative_r, adjacent_theta_idx + adjacent_relative_theta, z_idx + relative_z]
            end
        else
            # println("Calculated an invalid number of overlaps: $num_overlaps")
            # println("This is likely due to a cell somewhere being thinner than a particle")
            # Assign particle completely inside cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        end

        problem_remainder_volumes = findall(remainder_volumes .< -1e-11)
        
        # if !isempty(problem_remainder_volumes)
        #     println("WARNING: Total particle volume in cell calculated as negative in cell: ", cell_idx)
        #     println("Volume(s) calculated as: ", remainder_volumes[problem_remainder_volumes])
        #     println("$num_overlaps overlaps")
        #     case_counter[case + 1] += 1
        # end
        # total_case_counter[case + 1] += 1
        
        particle_volume_per_cell[cell_idx] += main_volume

        num_remainders = length(remainder_indices)
        for j in 1:num_remainders
            remainder_idx = remainder_indices[j]
            remainder_volume = remainder_volumes[j]
            # println(remainder_idx)
            remainder_cell_idx = find_global_cell_index_cylindrical(remainder_idx, r_divisions, z_divisions)
            particle_volume_per_cell[remainder_cell_idx] += remainder_volume
        end

    end
    # println(case_counter)
    # println(total_case_counter)
    return particle_volume_per_cell
end


function compute_centres_volume_per_cell_cylindrical(num_cells, num_particles, particle_volumes, r_data, theta_data, z_data, dr, r_divisions, dz, z_divisions, cylinder_base_level)
    # Initialise particle_volume_per_cell
    particle_volume_per_cell = zeros(Float64, num_cells)

    # Loop through all cells and calculate packing densities, then concentration of species 1
    for i in 1:num_particles
        r_idx, theta_idx, z_idx = compute_cell_index_cylindrical(r_data[i], theta_data[i], z_data[i], dr, r_divisions, dz, z_divisions, cylinder_base_level)
        if _is_index_outside_cylindrical(r_idx, z_idx, r_divisions, z_divisions)
            continue
        end
        cell_idx = find_global_cell_index_cylindrical([r_idx, theta_idx, z_idx], r_divisions, z_divisions)
        particle_volume_per_cell[cell_idx] += particle_volumes[i]
    end

    return particle_volume_per_cell
end



@inline function _particle_overlaps_cylindrical(
    cell_boundaries::Vector{Float64},
    r::Float64,
    theta::Float64,
    z::Float64,
    radius::Float64
)::Tuple{Vector{Float64}, BitVector}

    # Unpack boundaries
    r_min, r_max, theta_min, theta_max, z_min, z_max = cell_boundaries

    # Six signed distances: 
    # 1) below r_min, 2) above r_max
    # 3) “below” θ_min, 4) “above” θ_max  (converted to arc‐length)
    # 5) below z_min, 6) above z_max
    

    # Find which overlaps are valid
    if r_min < 0
        # Full cylindrical cell, don't check r_min or theta boundaries
        overlap_values = Float64[
            NaN,
            r     - r_max,
            NaN,
            NaN,
            z_min - z,
            z     - z_max
        ]

        overlaps = falses(6)
        overlaps[2] = abs(overlap_values[2]) < radius
        overlaps[5] = abs(overlap_values[5]) < radius
        overlaps[6] = abs(overlap_values[6]) < radius
    else
        overlap_values = Float64[
            r_min - r,
            r     - r_max,
            r * sin(angular_difference(theta,    theta_min)),
            r * sin(angular_difference(theta_max, theta   )),
            z_min - z,
            z     - z_max
        ]
        
        overlaps = abs.(overlap_values) .< radius
    end

    return overlap_values, overlaps
end


@inline function adjacent_radial_overlaps(r, theta, radius, r_idx, relative_r)
    adjacent_r_idx = r_idx + relative_r
    num_theta_divisions = 2*adjacent_r_idx - 1
    division_size = 2*pi/num_theta_divisions
    
    nearest_theta_boundary = round(theta / division_size) * division_size
    adjacent_theta_idx = Int(floor(theta / division_size)) + 1
    overlap_value = r*sin(angular_difference(nearest_theta_boundary, theta))
    
    if -radius < overlap_value < radius
        if overlap_value < 0
            # At a min boundary
            return overlap_value, adjacent_theta_idx, -1
        else
            # At a max boundary
            return overlap_value, adjacent_theta_idx, 1
        end
    else
        return nothing, adjacent_theta_idx, nothing
    end
end


"""
Guesses cylinder radius and base level assuming **at least one** particle is in 
contact with **bottom** and **edge** of the cylinder.
"""
function find_cylinder_parameters(; r_data=nothing,
                                    z_data=nothing,
                                    radii=nothing,
                                    mesh_vtk=nothing)
    if any(isnothing.([r_data, z_data, radii])) && isnothing(mesh_vtk)
        throw(ArgumentError("You stoopid. You didn't put any data in silly :("))
    end

    if length(r_data) != length(radii) || length(radii) != length(z_data)
        throw(ArgumentError("Arrays are different lengths stoopid"))
    end

    bottom_z_positions = z_data - radii
    cylinder_base_level = minimum(bottom_z_positions)

    outermost_r_positions = r_data + radii
    cylinder_radius = maximum(outermost_r_positions)

    return cylinder_radius, cylinder_base_level
end


@inline function compute_cell_index_cylindrical(r, theta, z, dr, r_divisions, dz, z_divisions, cylinder_base_level)
    # Compute z_idx
    z_idx = Int(floor((z - cylinder_base_level) / dz)) + 1
    z_idx = clamp(z_idx, 1, z_divisions)

    # Compute r_idx
    r_idx = Int(floor(r / dr)) + 1
    r_idx = clamp(r_idx, 1, r_divisions)

    # Compute theta_idx based on r_idx
    if r_idx == 1
        theta_idx = 1
    else
        # Number of theta divisions for this r_idx
        current_theta_div = 2 * r_idx - 1

        # Size of each theta bin
        d_theta = 2 * π / current_theta_div

        # Normalize theta to [0, 2π)
        theta_norm = mod(theta, 2 * π)

        # Compute theta_idx
        theta_idx = Int(floor(theta_norm / d_theta)) + 1
        theta_idx = clamp(theta_idx, 1, current_theta_div)
    end

    return r_idx, theta_idx, z_idx
end


@inline function find_global_cell_index_cylindrical(idx, r_divisions, z_divisions)
    r_idx, theta_idx, z_idx = idx
    # Compute cell_index within the current z layer
    if r_idx == 1
        cell_index_in_z = 1
    else
        # Cells are ordered with the inner cell first, followed by radial and angular cells
        # Sum of cells from inside to r_idx = r_idx^2
        # Therefore sum to r_idx - 1, then add theta_idx
        # theta_idx is modded to allow wraparounds
        cell_index_in_z = (r_idx - 1)^2 + mod(theta_idx - 1, 2*r_idx - 1) + 1
        if r_idx > r_divisions
            throw(ArgumentError("r_idx ($r_idx) greater than number of r_divisions ($r_divisions)"))
        end
        if z_idx > z_divisions
            throw(ArgumentError("z_idx ($z_idx) greater than number of z_divisions ($z_divisions)"))
        end
    end

    # Compute the global cell index
    global_cell_index = (z_idx - 1) * r_divisions^2 + cell_index_in_z

    return global_cell_index
end

@inline function _is_index_outside_cylindrical(r_idx, z_idx, r_divisions, z_divisions)
    return (r_idx > r_divisions) ||
           (z_idx < 1) || (z_idx > z_divisions)
end


end # module Cylindrical
