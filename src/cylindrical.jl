# cylindrical.jl
# This module contains all the necessary functions that are unique to
# calculating packing density in cylindrical coordinates

module Cylindrical

export compute_packing_cylindrical, calculate_segregation_intensity, calculate_lacey, calculate_lacey_old

# Import functions and structs from relevant modules
include("io.jl")
include("geometry.jl")
include("utils.jl")
include("cartesian.jl")
include("mesh.jl")

# Import specific functions from modules
using .IO: read_vtk_file, retrieve_coordinates

using .Geometry: convert_to_cylindrical,
                 calculate_angular_overlap_factor,
                 compute_cell_volume,
                 single_cap_intersection,
                 double_cap_intersection,
                 triple_cap_intersection,
                 sphere_cylinder_intersection,
                 sphere_cylinder_plane_intersection,
                 angular_difference

using .Utils: compute_automatic_boundaries,
              calculate_overlaps,
              calculate_active_overlap_values,
              is_inside_boundaries,
              is_outside_boundaries,
              centre_inside_boundaries,
              convert_boundaries_dictionary

using .Cartesian: calculate_particle_volume

using .MeshModule: Mesh,
                   get_mesh_boundaries,
                   get_total_cells,
                   get_cell_boundaries

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
function compute_packing_cylindrical(; file::Union{String, Nothing}=nothing,
                                      boundaries::Union{Vector{<:Real}, Dict{Symbol, <:Real}, Nothing}=nothing,
                                      r_data::Union{Vector{Float64}, Nothing}=nothing,
                                      theta_data::Union{Vector{Float64}, Nothing}=nothing,
                                      z_data::Union{Vector{Float64}, Nothing}=nothing,
                                      radii::Union{Vector{Float64}, Nothing}=nothing,
                                      accurate_cylindrical::Union{Bool, Nothing}=nothing)::Float64

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

    r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

    if r_min == 0.0 && abs(theta_max - theta_min - 2π) < 1e-8
        throw(ArgumentError("r_min = 0 with a full circular theta range is invalid. Use r_min < 0 for a full cylinder."))
    end

    full_cylindrical_cell = r_min < 0.0

    # Step 3: Calculate angular overlap factor
    factor = calculate_angular_overlap_factor(r_data, radii)

    # Step 4: Calculate overlaps and active overlap values
    overlaps = calculate_overlaps(; r_data=r_data, theta_data=theta_data, z_data=z_data,
                                   radii=radii, boundaries=boundaries, factor=factor, system=:cylindrical)

    total_particles = length(radii)
    active_overlap_values = calculate_active_overlap_values(total_particles;
                                                            x_data=nothing,
                                                            y_data=nothing,
                                                            r_data=r_data,
                                                            theta_data=theta_data,
                                                            z_data=z_data,
                                                            boundaries=boundaries,
                                                            overlaps=overlaps,
                                                            system=:cylindrical)

    # Step 5: Precompute particle volumes
    full_particle_volumes = (4 / 3) * π .* (radii .^ 3)

    # Step 6: Determine masks for particles completely inside or outside boundaries
    inside_mask = is_inside_boundaries(; x_data=nothing, y_data=nothing,
                                         r_data=r_data, theta_data=theta_data, z_data=z_data,
                                         boundaries=boundaries, radii=radii,
                                         factor=factor,
                                         system=:cylindrical)

    outside_mask = is_outside_boundaries(; x_data=nothing, y_data=nothing,
                                           r_data=r_data, theta_data=theta_data, z_data=z_data,
                                           boundaries=boundaries, radii=radii,
                                           factor=factor,
                                           system=:cylindrical)

    # Enable accurate cylindrical calculation if the full cell condition is met
    if isnothing(accurate_cylindrical) && full_cylindrical_cell
        accurate_cylindrical = true
    end

    # Step 7: Initialize total particle volume
    total_particle_volume = sum(full_particle_volumes[inside_mask])

    # Skip fully outside particles and process only "neither" cases
    neither_mask = .~(inside_mask .| outside_mask)
    indices_neither = findall(neither_mask)

    if !isempty(indices_neither)
        if accurate_cylindrical
            # Calculate partial volumes accurately for overlapping particles
            total_particle_volume += sum(calculate_partial_volumes_cyl(i, radii, active_overlap_values, r_data)
                                         for i in indices_neither)
        else
            # Calculate partial volumes using planar approximation for overlapping particles
            total_particle_volume += sum(calculate_particle_volume(i, radii, active_overlap_values)
                                         for i in indices_neither)
        end
    end

    # Step 8: Compute cell volume
    cell_volume = compute_cell_volume(; boundaries=boundaries, system=:cylindrical)

    # Step 9: Calculate packing density
    packing_density = total_particle_volume / cell_volume

    return packing_density

end


"""
    calculate_particle_volume_cyl(i::Int, radii::Vector{Float64}, 
                                  active_overlap_values::Matrix{Float64}, 
                                  r_data::Vector{Float64})::Float64

Calculate the volume contribution of a particle in a cylindrical system.

# Arguments
- `i::Int`: Index of the particle being evaluated.
- `radii::Vector{Float64}`: Radii of the particles.
- `active_overlap_values::Matrix{Float64}`: Overlap distances with boundaries.
- `r_data::Vector{Float64}`: Radial positions of the particles.

# Returns
- `Float64`: The volume contribution of the particle.
"""
function calculate_partial_volumes_cyl(i::Int,
                                       radii::Vector{Float64},
                                       active_overlap_values::Matrix{Float64},
                                       r_data::Vector{Float64}
                                       )::Float64
    # Initialize partial volume
    partial_volume = 0.0

    # Extract overlaps using nested ternary operators
    r_overlap = !isnan(active_overlap_values[i, 1]) ? active_overlap_values[i, 1] :
                (!isnan(active_overlap_values[i, 2]) ? active_overlap_values[i, 2] : nothing)

    theta_overlap = !isnan(active_overlap_values[i, 3]) ? active_overlap_values[i, 3] :
                   (!isnan(active_overlap_values[i, 4]) ? active_overlap_values[i, 4] : nothing)

    z_overlap = !isnan(active_overlap_values[i, 5]) ? active_overlap_values[i, 5] :
                (!isnan(active_overlap_values[i, 6]) ? active_overlap_values[i, 6] : nothing)

    # Determine radial overlap type and boundary
    overlaps_r, min_boundary_r = if !isnan(active_overlap_values[i, 1])
        (true, true)  # Overlaps with r_min
    elseif !isnan(active_overlap_values[i, 2])
        (true, false)  # Overlaps with r_max
    else
        (false, false)
    end

    # Check for overlaps in theta and z dimensions
    overlaps_theta = !isnothing(theta_overlap)
    overlaps_z = !isnothing(z_overlap)

    # Count the number of overlaps
    num_overlaps = sum([overlaps_r, overlaps_theta, overlaps_z])

    # Calculate partial volume based on the number of overlaps
    if num_overlaps == 1
        if overlaps_r
            partial_volume = sphere_cylinder_intersection(
                radii[i], r_data[i], r_overlap; min_boundary=min_boundary_r
            )
        elseif overlaps_theta
            partial_volume = single_cap_intersection(radii[i], theta_overlap)
        elseif overlaps_z
            partial_volume = single_cap_intersection(radii[i], z_overlap)
        end
    elseif num_overlaps == 2
        if overlaps_r && overlaps_theta
            partial_volume = double_cap_intersection(radii[i], r_overlap, theta_overlap)
        elseif overlaps_r && overlaps_z
            partial_volume = sphere_cylinder_plane_intersection(
                radii[i], r_data[i], r_overlap, z_overlap; min_boundary=min_boundary_r
            )
        elseif overlaps_theta && overlaps_z
            partial_volume = double_cap_intersection(radii[i], theta_overlap, z_overlap)
        end
    elseif num_overlaps == 3
        if overlaps_r && overlaps_theta && overlaps_z
            partial_volume = triple_cap_intersection(
                radii[i], r_overlap, theta_overlap, z_overlap
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
function calculate_segregation_intensity(
    data_1::Dict, data_2::Dict, 
    cylinder_radius::Float64,
    cylinder_base_level::Float64, 
    cylinder_height::Float64;
    target_num_cells::Int=200,
    packing_threshold::Float64=0.05,
    calculate_partial_volumes::Bool=true
)::Float64
    # Step 1: Determine mesh divisions
    z_divisions = max(1, round(Int, target_num_cells^(1 / 3)))
    num_cells_slice = target_num_cells / z_divisions
    theta_divisions = 3  # Fixed as per original code
    r_divisions = max(1, round(Int, sqrt(num_cells_slice)))

    # Step 2: Extract particle data from both datasets
    x_data_1, y_data_1, z_data_1, radii_1 = retrieve_coordinates(data_1)
    x_data_2, y_data_2, z_data_2, radii_2 = retrieve_coordinates(data_2)

    # Convert Cartesian to cylindrical coordinates
    r_data_1, theta_data_1 = convert_to_cylindrical(x_data_1, y_data_1)
    r_data_2, theta_data_2 = convert_to_cylindrical(x_data_2, y_data_2)

    factor_1 = calculate_angular_overlap_factor(r_data_1, radii_1)
    factor_2 = calculate_angular_overlap_factor(r_data_2, radii_2)

    # Step 3: Compute total volumes of both datasets
    particle_volumes_1 = (4 / 3) * π * (radii_1 .^ 3)
    particle_volumes_2 = (4 / 3) * π * (radii_2 .^ 3)
    total_volume_1 = sum(particle_volumes_1)
    total_volume_2 = sum(particle_volumes_2)
    conc_mean = total_volume_1 / (total_volume_1 + total_volume_2)

    # Step 4: Create the cylindrical mesh
    divisions = Dict("r" => r_divisions, "theta" => theta_divisions, "z" => z_divisions)
    params = Dict("cylinder_radius" => cylinder_radius, 
                  "cylinder_base_level" => cylinder_base_level, 
                  "cylinder_height" => cylinder_height)

    cylindrical_mesh = Mesh(:cylindrical, divisions; params=params)
    mesh_boundaries = get_mesh_boundaries(cylindrical_mesh)
    num_cells = get_total_cells(cylindrical_mesh)

    # Step 5: Initialize sparse arrays for packing densities
    concs_1_sparse = fill(NaN, num_cells)

    # Compute cell volume using the first cell's boundaries
    cell_volume = pi * (cylinder_radius / r_divisions)^2 * (cylinder_height / z_divisions)


    # Step 6: Decide whether to calculate partial volumes or not
    if calculate_partial_volumes
        # Loop through all cells and calculate packing densities, then concentration of species 1
        for index in 1:num_cells
            boundaries = mesh_boundaries[index, :]

            # Compute packing densities for data_1 and data_2
            packing_density_1 = compute_packing_cylindrical(; boundaries=boundaries, 
                                                              r_data=r_data_1, 
                                                              theta_data=theta_data_1, 
                                                              z_data=z_data_1, 
                                                              radii=radii_1, 
                                                              accurate_cylindrical=false)

            packing_density_2 = compute_packing_cylindrical(; boundaries=boundaries, 
                                                              r_data=r_data_2, 
                                                              theta_data=theta_data_2, 
                                                              z_data=z_data_2, 
                                                              radii=radii_2, 
                                                              accurate_cylindrical=false)

            # Total packing density
            packing_density_total = packing_density_1 + packing_density_2

            # Exclude cells with low packing density
            if packing_density_total > packing_threshold
                concs_1_sparse[index] = packing_density_1 / packing_density_total
            else
                concs_1_sparse[index] = NaN
            end
        end
    else
        # Binning and Aggregation

        # Step 6.1: Define bin sizes
        dr = cylinder_radius / r_divisions
        dz = cylinder_height / z_divisions

        # Step 6.2: Helper Function to Compute Linear Cell Index
        function compute_cell_index(r, theta, z)
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

            # Compute cell_index within the current z layer
            if r_idx == 1
                cell_index_in_z = 1
            else
                # Cells are ordered with the inner cell first, followed by radial and angular cells
                # cell_index_in_z = 1 + sum_{i=2}^{r_idx} (2*i -1) + theta_idx
                # Since sum_{i=2}^{r_idx} (2*i -1) = r_idx^2 -1
                cell_index_in_z = r_idx^2 + theta_idx
            end

            # Compute the global cell index
            global_cell_index = (z_idx - 1) * r_divisions^2 + cell_index_in_z

            return global_cell_index
        end

        # Step 6.3: Assign Particles to Cells for Both Datasets
        # Preallocate arrays for cell indices
        num_particles_1 = length(r_data_1)
        cell_indices_1 = Vector{Int}(undef, num_particles_1)

        num_particles_2 = length(r_data_2)
        cell_indices_2 = Vector{Int}(undef, num_particles_2)

        # Assign cell indices for data_1
        for i in 1:num_particles_1
            cell_indices_1[i] = compute_cell_index(r_data_1[i], theta_data_1[i], z_data_1[i])
        end

        # Assign cell indices for data_2
        for i in 1:num_particles_2
            cell_indices_2[i] = compute_cell_index(r_data_2[i], theta_data_2[i], z_data_2[i])
        end

        # Step 6.4: Aggregate Volumes per Cell
        volume_per_cell_1 = zeros(Float64, num_cells)
        volume_per_cell_2 = zeros(Float64, num_cells)

        for i in 1:length(cell_indices_1)
            idx = cell_indices_1[i]
            # Ensure the cell index is within valid range
            if idx >=1 && idx <= num_cells
                volume_per_cell_1[idx] += particle_volumes_1[i]
            else
                # Optionally, handle particles outside the mesh boundaries
                # For now, ignore them
            end
        end

        for i in 1:length(cell_indices_2)
            idx = cell_indices_2[i]
            if idx >=1 && idx <= num_cells
                volume_per_cell_2[idx] += particle_volumes_2[i]
            else
                # Optionally, handle particles outside the mesh boundaries
                # For now, ignore them
            end
        end

        # Step 6.5: Compute Total Packing Density per Cell
        total_particle_volume_per_cell = (volume_per_cell_1 .+ volume_per_cell_2)

        # Step 6.6: Compute Concentration of Species 1
        # Initialize with NaN
        concs_1_sparse = fill(NaN, num_cells)

        # Determine which cells meet the packing density threshold
        valid_cells_mask = total_particle_volume_per_cell .> packing_threshold .* cell_volume

        # Calculate concentration for valid cells
        concs_1_sparse[valid_cells_mask] = volume_per_cell_1[valid_cells_mask] ./ total_particle_volume_per_cell[valid_cells_mask]

    end

    # Step 7: Remove NaNs and count valid cells
    concs_1_filtered = filter(!isnan, concs_1_sparse)
    num_valid_cells = length(concs_1_filtered)

    num_invalid_cells = num_cells - num_valid_cells

    if num_invalid_cells > 0.5 * num_cells
        println("Unusually high number of invalid cells: $(num_invalid_cells) out of $(num_cells)")
    end

    if num_valid_cells == 0
        return NaN  # No valid cells
    end

    # Step 8: Calculate segregation intensity
    numerator = sum((concs_1_filtered .- conc_mean).^2)
    I_S_max = sqrt(conc_mean * (1 - conc_mean))
    segregation_intensity = sqrt(numerator / num_valid_cells) / I_S_max

    return segregation_intensity
end


function calculate_lacey_old(data_1::Dict,
                             data_2::Dict,
                             cylinder_radius::Float64,
                             cylinder_base_level::Float64,
                             cylinder_height::Float64,
                             ;
                             target_num_cells::Int=200,
                             packing_threshold::Float64=0.05,
                             output_num_cells::Bool=false,
                             calculate_partial_volumes::Bool=true,
                             clamp_0_to_1::Bool=true)::Float64
    # Step 1: Determine mesh divisions
    z_divisions = Int(round(target_num_cells^(1 / 3)))
    num_cells_slice = target_num_cells / z_divisions
    theta_divisions = 3
    r_divisions = Int(round(sqrt(num_cells_slice)))

    num_cells = Int(round(z_divisions * r_divisions^2))
    if output_num_cells
        println("Target number of cells: $target_num_cells, Actual number: $num_cells")
    end

    radius_inner = cylinder_radius / r_divisions
    cell_volume = π * radius_inner^2 * cylinder_height / z_divisions

    # Step 2: Extract particle data from both datasets
    x_data_1, y_data_1, z_data_1, radii_1 = retrieve_coordinates(data_1)
    x_data_2, y_data_2, z_data_2, radii_2 = retrieve_coordinates(data_2)

    # Convert Cartesian to cylindrical coordinates
    r_data_1, theta_data_1 = convert_to_cylindrical(x_data_1, y_data_1)
    r_data_2, theta_data_2 = convert_to_cylindrical(x_data_2, y_data_2)

    # Step 3: Create the cylindrical mesh
    divisions = Dict("r" => r_divisions, "theta" => theta_divisions, "z" => z_divisions)
    params = Dict("cylinder_radius" => cylinder_radius,
                  "cylinder_base_level" => cylinder_base_level,
                  "cylinder_height" => cylinder_height)

    cylindrical_mesh = Mesh(:cylindrical, divisions; params=params)
    mesh_boundaries = cylindrical_mesh.cell_boundaries
    num_cells = cylindrical_mesh.total_cells

    # Step 4: Initialize sparse arrays
    concs_1_sparse = fill(NaN, num_cells)
    particle_volumes_sparse = fill(NaN, num_cells)

    # Step 5: Compute packing densities and volumes for each cell
    for index in 1:num_cells
        boundaries = mesh_boundaries[index, :]

        # Compute packing densities for data_1 and data_2
        packing_density_1 = compute_packing_cylindrical(; boundaries=boundaries,
                                                        r_data=r_data_1,
                                                        theta_data=theta_data_1,
                                                        z_data=z_data_1,
                                                        radii=radii_1,
                                                        accurate_cylindrical=false)

        packing_density_2 = compute_packing_cylindrical(; boundaries=boundaries,
                                                        r_data=r_data_2,
                                                        theta_data=theta_data_2,
                                                        z_data=z_data_2,
                                                        radii=radii_2,
                                                        accurate_cylindrical=false)

        packing_density_total = packing_density_1 + packing_density_2

        # Exclude cells with low packing density
        if packing_density_total > packing_threshold
            conc_1 = packing_density_1 / packing_density_total
            concs_1_sparse[index] = conc_1
            particle_volumes_sparse[index] = packing_density_total * cell_volume
        end
    end

    # Step 6: Filter valid cells
    concs_1 = filter(!isnan, concs_1_sparse)
    particle_volumes = filter(!isnan, particle_volumes_sparse)

    num_valid_cells = length(concs_1)

    # Step 7: Compute effective particles per cell
    eff_particles_per_cell = sum(particle_volumes)^2 / sum(particle_volumes.^2)

    # Step 8: Compute weighted concentration mean
    weighted_conc_mean = sum(concs_1 .* particle_volumes) / sum(particle_volumes)

    # Step 9: Handle edge cases
    if num_valid_cells <= 1
        return NaN
    end

    # Step 10: Compute Lacey index
    # Actual variance of concentrations
    S_actual = sum((concs_1 .- weighted_conc_mean).^2) / num_valid_cells

    # Random variance (binomial distribution assumption)
    S_random = weighted_conc_mean * (1 - weighted_conc_mean) / eff_particles_per_cell

    # Maximum variance (complete segregation)
    S_maximum = weighted_conc_mean * (1 - weighted_conc_mean)

    # Final calculation of Lacey index
    M = (S_maximum - S_actual) / (S_maximum - S_random)

    return M
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
function calculate_lacey(data_1::Dict,
                         data_2::Dict
                         ;
                         cylinder_radius::Union{<:Real, Nothing}=nothing,
                         cylinder_base_level::Union{<:Real, Nothing}=nothing,
                         cylinder_height::Union{<:Real, Nothing}=nothing,
                         target_num_cells::Int=200,
                         packing_threshold::Real=0.05,
                         output_num_cells::Bool=false,
                         calculate_partial_volumes::Bool=true,
                         clamp_0_to_1::Bool=true)::Float64

    # Step 1: Determine mesh divisions
    z_divisions = max(1, round(Int, target_num_cells^(1 / 3)))
    num_cells_slice = target_num_cells / z_divisions
    theta_divisions = 3  # Fixed as per original code
    r_divisions = max(1, round(Int, sqrt(num_cells_slice)))

    # Step 2: Extract particle data from both datasets
    x_data_1, y_data_1, z_data_1, radii_1 = retrieve_coordinates(data_1)
    x_data_2, y_data_2, z_data_2, radii_2 = retrieve_coordinates(data_2)

    # Convert Cartesian to cylindrical coordinates
    r_data_1, theta_data_1 = convert_to_cylindrical(x_data_1, y_data_1)
    r_data_2, theta_data_2 = convert_to_cylindrical(x_data_2, y_data_2)

    if isnothing(cylinder_height)
        throw(ArgumentError("Cylinder height must be provided stoopid"))
    end

    estimated_cylinder_radius, estimated_cylinder_base_level = (
        find_cylinder_parameters(; r_data=vcat(r_data_1, r_data_2),
                                   z_data=vcat(z_data_1, z_data_2),
                                   radii=vcat(radii_1, radii_2))
    )

    if isnothing(cylinder_radius)
        cylinder_radius = estimated_cylinder_radius
        # println("No cylinder radius provided, using esimation: $(cylinder_radius)")
    end

    if isnothing(cylinder_base_level)
        cylinder_base_level = estimated_cylinder_base_level
        # println("No cylinder base level provided, using esimation: $(cylinder_base_level)")
    end

    # Step 3: Compute total volumes of both datasets
    particle_volumes_1 = (4 / 3) * π * (radii_1 .^ 3)
    particle_volumes_2 = (4 / 3) * π * (radii_2 .^ 3)

    # Step 4: Create the cylindrical mesh
    divisions = Dict("r" => r_divisions, "theta" => theta_divisions, "z" => z_divisions)
    params = Dict("cylinder_radius" => cylinder_radius, 
                  "cylinder_base_level" => cylinder_base_level, 
                  "cylinder_height" => cylinder_height)

    cylindrical_mesh = Mesh(:cylindrical, divisions; params=params)
    mesh_boundaries = get_mesh_boundaries(cylindrical_mesh)
    num_cells = get_total_cells(cylindrical_mesh)
    # num_cells = r_divisions^2 * z_divisions

    num_particles_1 = length(r_data_1)
    num_particles_2 = length(r_data_2)
    num_particles = 

    volume_per_cell_1 = zeros(Float64, num_cells)
    volume_per_cell_2 = zeros(Float64, num_cells)

    if output_num_cells
        println("Target number of cells: $target_num_cells")
        println("Actual number of cells: $num_cells")
        println("Radial divisions:       $r_divisions")
        println("Axial divisions:        $z_divisions")
    end

    # Step 5: Initialize sparse arrays for packing densities
    concs_1_sparse = fill(NaN, num_cells)
    total_particle_volume_per_cell = fill(NaN, num_cells)

    # Compute cell volume using the first cell's boundaries
    cell_volume = pi * (cylinder_radius / r_divisions)^2 * (cylinder_height / z_divisions)

    dr = cylinder_radius / r_divisions
    dz = cylinder_height / z_divisions

    # Step 6: Decide whether to calculate partial volumes or not
    if calculate_partial_volumes
        volume_per_cell_1 = compute_volume_per_cell(mesh_boundaries, num_cells, num_particles_1, particle_volumes_1, r_data_1, theta_data_1, z_data_1, radii_1, dr, r_divisions, dz, z_divisions, cylinder_base_level)
        volume_per_cell_2 = compute_volume_per_cell(mesh_boundaries, num_cells, num_particles_2, particle_volumes_2, r_data_2, theta_data_2, z_data_2, radii_2, dr, r_divisions, dz, z_divisions, cylinder_base_level)
    else
        # Binning and Aggregation

        # Step 6.3: Assign Particles to Cells for Both Datasets

        # Assign cell indices for data_1
        for i in 1:num_particles_1
            r_idx, theta_idx, z_idx = compute_cell_index(r_data_1[i], theta_data_1[i], z_data_1[i], dr, r_divisions, dz, z_divisions, cylinder_base_level)
            cell_idx = find_global_cell_index([r_idx, theta_idx, z_idx], r_divisions, z_divisions)
            volume_per_cell_1[cell_idx] += particle_volumes_1[i]
        end

        # Assign cell indices for data_2
        for i in 1:num_particles_2
            r_idx, theta_idx, z_idx = compute_cell_index(r_data_2[i], theta_data_2[i], z_data_2[i], dr, r_divisions, dz, z_divisions, cylinder_base_level)
            cell_idx = find_global_cell_index([r_idx, theta_idx, z_idx], r_divisions, z_divisions)
            volume_per_cell_2[cell_idx] += particle_volumes_2[i]
        end
    end

    # Step 6.5: Compute Total Packing Density per Cell
    total_particle_volume_per_cell = (volume_per_cell_1 .+ volume_per_cell_2)
    
    # println("Total of assigned cell volumes: $(sum(total_particle_volume_per_cell))")
    # println("Actual total particle volume:   $(sum(particle_volumes_1) + sum(particle_volumes_2))")

    conc_mean_valid = sum(volume_per_cell_1) / sum(total_particle_volume_per_cell)

    # Step 6.6: Compute Concentration of Species 1
    # Initialize with NaN
    concs_1_sparse = fill(NaN, num_cells)

    # Determine which cells meet the packing density threshold
    valid_cells_mask = total_particle_volume_per_cell .> packing_threshold .* cell_volume

    # Calculate concentration for valid cells
    concs_1_sparse[valid_cells_mask] = volume_per_cell_1[valid_cells_mask] ./ total_particle_volume_per_cell[valid_cells_mask]
    
    particle_volumes = total_particle_volume_per_cell[valid_cells_mask]

    # Step 6: Filter valid cells
    concs_1_valid = filter(!isnan, concs_1_sparse)
    
    num_valid_cells = length(concs_1_valid)

    # Step 7: Compute effective particles per cell
    eff_particles_per_cell = sum(particle_volumes)^2 / sum(particle_volumes.^2)

    # Step 8: Compute weighted concentration mean
    weighted_conc_mean = sum(concs_1_valid .* particle_volumes) / sum(particle_volumes)

    # Step 9: Handle edge cases
    if num_valid_cells <= 1
        return NaN
    end

    # Step 10: Compute Lacey index
    # Actual variance of concentrations
    S_actual = sum((concs_1_valid .- weighted_conc_mean).^2) / num_valid_cells

    # Random variance (binomial distribution assumption)
    S_random = weighted_conc_mean * (1 - weighted_conc_mean) / eff_particles_per_cell
    
    # Maximum variance (complete segregation, binomial distribution assumption)
    S_maximum = weighted_conc_mean * (1 - weighted_conc_mean)

    # Actual variance of concentrations
    # S_actual = sum((concs_1_valid .- conc_mean_valid).^2) / num_valid_cells

    # # Random variance: Fully mixed -> All cells have a concentration of conc_mean, therefore no variance
    # S_random = conc_mean_valid * (1 - conc_mean_valid) / num_valid_cells

    # # Maximum variance (complete segregation, binomial distribution assumption)
    # S_maximum = conc_mean_valid * (1 - conc_mean_valid)

    # println("S_random: $S_random")
    # println("S_maximum: $S_maximum")
    # println("S_actual: $S_actual")


    # Final calculation of Lacey index
    M = (S_maximum - S_actual) / (S_maximum - S_random)

    if clamp_0_to_1
        M = clamp(M, 0, 1)
    end

    return M
end


function compute_volume_per_cell(mesh_boundaries, num_cells, num_particles, particle_volumes, r_data, theta_data, z_data, radii, dr, r_divisions, dz, z_divisions, cylinder_base_level)
    # Initialise particle_volume_per_cell
    particle_volume_per_cell = zeros(Float64, num_cells)
    
    # Loop through all cells and calculate packing densities, then concentration of species 1
    for i in 1:num_particles
        r_idx, theta_idx, z_idx = compute_cell_index(r_data[i], theta_data[i], z_data[i], dr, r_divisions, dz, z_divisions, cylinder_base_level)
        cell_idx = find_global_cell_index([r_idx, theta_idx, z_idx], r_divisions, z_divisions)
        overlap_values, overlaps = particle_overlaps(mesh_boundaries[cell_idx, :], r_data[i], theta_data[i], z_data[i], radii[i])
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

        if num_overlaps == 0
            # Particle completely inside cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        elseif num_overlaps == 1
            overlap_value = overlap_values[overlaps][1]
            if overlaps_r
                main_volume = sphere_cylinder_intersection(radii[i], r_data[i], overlap_value; min_boundary=min_boundary)
                if isnothing(adjacent_overlap)
                    remainder_volumes = [particle_volumes[i] - main_volume]
                    remainder_indices = [[r_idx + relative_r, adjacent_theta_idx, z_idx]]
                else
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
                    remainder_indices = fill(zeros(Int, 2), 2)
                    remainder_volumes = zeros(Float64, 2)
                    remainder_volumes[1] = sphere_cylinder_intersection(radii[i], r_data[i], -overlap_vals[1]; min_boundary=!min_boundary)
                    remainder_indices[1] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                    remainder_volumes[2] = particle_volumes[i] - main_volume - remainder_volumes[1]
                    remainder_indices[2] = [r_idx, theta_idx + relative_theta, z_idx]
                else
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
                    remainder_indices = fill(zeros(Int, 3), 3)
                    remainder_volumes = zeros(Float64, 3)
                    remainder_volumes[1] = sphere_cylinder_intersection(radii[i], r_data[i], overlap_vals[1]) - main_volume
                    remainder_indices[1] = [r_idx, theta_idx, z_idx + relative_z]
                    remainder_volumes[2] = double_cap_intersection(radii[i], -overlap_vals[1], overlap_vals[2])
                    remainder_indices[2] = [r_idx + relative_r, adjacent_theta_idx, z_idx]
                    remainder_volumes[3] = particle_volumes[i] - main_volume - remainder_volumes[1] - remainder_volumes[2]
                    remainder_indices[3] = [r_idx + relative_r, adjacent_theta_idx, z_idx + relative_z]
                else
                    remainder_indices = fill(zeros(Int, 5), 5)
                    remainder_volumes = zeros(Float64, 5)
                    remainder_volumes[1] = sphere_cylinder_intersection(radii[i], r_data[i], overlap_vals[1]) - main_volume
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
                println("Calculated 2 overlaps but could not match the case")
                println("Overlap values: $overlap_values")
                println("Overlaps mask: $overlaps")
                println("Mesh boundaries: $(mesh_boundaries[cell_idx, :])")
                println("Coordinates: $(r_data[i]), $(theta_data[i]), $(z_data[i])")
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
            println("Calculated an invalid number of overlaps: $num_overlaps")
        end
        

        
        particle_volume_per_cell[cell_idx] += main_volume

        num_remainders = length(remainder_indices)
        for j in 1:num_remainders
            remainder_idx = remainder_indices[j]
            remainder_volume = remainder_volumes[j]
            # println(remainder_idx)
            remainder_cell_idx = find_global_cell_index(remainder_idx, r_divisions, z_divisions)
            particle_volume_per_cell[remainder_cell_idx] += remainder_volume
        end

    end

    return particle_volume_per_cell
end


@inline function particle_overlaps(cell_boundaries, r, theta, z, radius)
    r_min, r_max, theta_min, theta_max, z_min, z_max = cell_boundaries
    tolerance = 1e-10
    overlap_values = [
        r_min - r,
        r - r_max,
        r*sin(angular_difference(theta, theta_min)),
        r*sin(angular_difference(theta_max, theta)),
        z_min - z,
        z - z_max
    ]
    if r_min < 0
        overlaps = [
            false,
            -radius - tolerance < overlap_values[2] < radius + tolerance,
            false,
            false,
            -radius - tolerance < overlap_values[5] < radius + tolerance,
            -radius - tolerance < overlap_values[6] < radius + tolerance,
        ]
    else
        overlaps = -radius - tolerance .< overlap_values .< radius + tolerance
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
# I love my wife #
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


@inline function compute_cell_index(r, theta, z, dr, r_divisions, dz, z_divisions, cylinder_base_level)
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


@inline function find_global_cell_index(idx, r_divisions, z_divisions)
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


end # module Cylindrical
