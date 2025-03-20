# cartesian.jl
# This module contains all the necessary functions that are unique to
# calculating packing density in cartesian coordinates

module Cartesian

# Import functions and structs from relevant modules
# include("io.jl")
# include("geometry.jl")
# include("utils.jl")
# include("mesh.jl")

# Using the modules (if functions are exported)
# using .IO, .Geometry, .Utils, .Cartesian, .MeshModule

# Import specific functions from modules
using ..IO: read_vtk_file, retrieve_coordinates

using ..Geometry: compute_cell_volume,
                 single_cap_intersection,
                 double_cap_intersection,
                 triple_cap_intersection

using ..Utils: compute_automatic_boundaries,
              calculate_overlaps,
              calculate_active_overlap_values,
              is_inside_boundaries,
              is_outside_boundaries,
              convert_boundaries_dictionary

using ..MeshModule: Mesh,
                   get_mesh_boundaries,
                   get_total_cells,
                   get_cell_boundaries,
                   compute_divisions


export _compute_packing_cartesian, _calculate_lacey_cartesian, _calculate_segregation_intensity_cartesian, _compute_volume_per_cell_cartesian


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
function _compute_packing_cartesian(; file::Union{String, Nothing}=nothing,
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


function _calculate_segregation_intensity_cartesian(
        data_1::Dict,
        data_2::Dict
        ;
        mesh::Union{Mesh, Nothing}=nothing,
        calculate_partial_volumes::Bool=true,
        clamp_0_to_1::Bool=true,
        verbose::Bool=false
    )::Float64

    volume_per_cell_1, volume_per_cell_2, real_pv_1, real_pv_2 = _compute_volume_per_cell_cartesian(
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


function _calculate_lacey_cartesian(
        data_1::Dict,
        data_2::Dict
        ;
        mesh::Union{Mesh, Nothing}=nothing,
        calculate_partial_volumes::Union{Bool, Nothing}=true,
        clamp_0_to_1::Union{Bool, Nothing}=true,
        verbose::Union{Bool, Nothing}=false
    )::Float64

    volume_per_cell_1, volume_per_cell_2, real_pv_1, real_pv_2 = _compute_volume_per_cell_cartesian(
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


function _compute_volume_per_cell_cartesian(data_1::Dict,
                                 data_2::Dict
                                 ;
                                 mesh::Union{Mesh, Nothing}=nothing,
                                 calculate_partial_volumes::Union{Bool, Nothing}=true,
                                 verbose::Union{Bool, Nothing}=false)

    # Extract mesh information
    params = mesh.params
    p = (; params...)
    divisions = mesh.divisions
    x_min, x_max, y_min, y_max, z_min, z_max = p.x_min, p.x_max, p.y_min, p.y_max, p.z_min, p.z_max
    L_x, L_y, L_z = (
        x_max - x_min,
        y_max - y_min,
        z_max - z_min
    )

    # Extract particle data from both datasets
    x_data_1, y_data_1, z_data_1, radii_1 = retrieve_coordinates(data_1)
    x_data_2, y_data_2, z_data_2, radii_2 = retrieve_coordinates(data_2)

    # Step 3: Compute total volumes of both datasets
    particle_volumes_1 = (4 / 3) * π * (radii_1 .^ 3)
    particle_volumes_2 = (4 / 3) * π * (radii_2 .^ 3)

    max_particle_volume = maximum(vcat(particle_volumes_1, particle_volumes_2))

    # Calculate division sizes and individual cell volume
    division_vals = collect(values(divisions))
    dx, dy, dz = [L_x, L_y, L_z] ./ division_vals
    cell_volume = dx * dy * dz
    recip_dx, recip_dy, recip_dz = 1 ./ [dx, dy, dz]

    if calculate_partial_volumes && (cell_volume < max_particle_volume)
        if verbose println("WARNING: Mesh resolution is too fine; refusing to calculate partial volumes") end
        calculate_partial_volumes = false
    end

    # Step 4: Create the cartesian mesh
    mesh_boundaries = get_mesh_boundaries(mesh)
    num_cells = get_total_cells(mesh)

    num_particles_1 = length(x_data_1)
    num_particles_2 = length(x_data_2)

    volume_per_cell_1 = zeros(Float64, num_cells)
    volume_per_cell_2 = zeros(Float64, num_cells)

    # Step 5: Initialise zeros array for cell volumes
    total_particle_volume_per_cell = fill(0.0, num_cells)

    # Step 6: Decide whether to calculate partial volumes or not
    if calculate_partial_volumes
        # Calculate partial volume contribution to overlapped cells
        # Refer to partial volume calculation function
        volume_per_cell_1 = _compute_partial_volume_per_cell_cartesian(mesh_boundaries, particle_volumes_1, num_cells, num_particles_1, x_data_1, y_data_1, z_data_1, radii_1, x_min, y_min, z_min, recip_dx, recip_dy, recip_dz, division_vals[1], division_vals[2], division_vals[3])
        volume_per_cell_2 = _compute_partial_volume_per_cell_cartesian(mesh_boundaries, particle_volumes_2, num_cells, num_particles_2, x_data_2, y_data_2, z_data_2, radii_2, x_min, y_min, z_min, recip_dx, recip_dy, recip_dz, division_vals[1], division_vals[2], division_vals[3])
    else
        # Not calculating partial volumes
        # Bin particles based on their centre positions

        # Assign cell indices for data_1
        for i in 1:num_particles_1
            if _is_particle_outside(x_data_1[i], y_data_1[i], z_data_1[i], p)
                # println("Particle found outside specified boundaries")
                continue
            end
            x_idx, y_idx, z_idx = _compute_cell_index_cartesian(x_data_1[i], y_data_1[i], z_data_1[i], x_min, y_min, z_min, recip_dx, recip_dy, recip_dz)
            cell_idx = _find_global_cell_index_cartesian(x_idx, y_idx, z_idx, division_vals[1], division_vals[2])
            volume_per_cell_1[cell_idx] += particle_volumes_1[i]
        end

        # Assign cell indices for data_2
        for i in 1:num_particles_2
            if _is_particle_outside(x_data_2[i], y_data_2[i], z_data_2[i], p)
                # println("Particle found outside specified boundaries")
                continue
            end
            x_idx, y_idx, z_idx = _compute_cell_index_cartesian(x_data_2[i], y_data_2[i], z_data_2[i], x_min, y_min, z_min, recip_dx, recip_dy, recip_dz)
            cell_idx = _find_global_cell_index_cartesian(x_idx, y_idx, z_idx, division_vals[1], division_vals[2])
            volume_per_cell_2[cell_idx] += particle_volumes_2[i]
        end
    end

    return volume_per_cell_1, volume_per_cell_2, particle_volumes_1, particle_volumes_2
end


function _compute_partial_volume_per_cell_cartesian(mesh_boundaries, particle_volumes, num_cells, num_particles, x_data, y_data, z_data, radii, x_min, y_min, z_min, recip_dx, recip_dy, recip_dz, x_divisions, y_divisions, z_divisions)
    # Initialise particle_volume_per_cell
    particle_volume_per_cell = zeros(Float64, num_cells)

    # case_counter = zeros(Int, 9)
    # total_case_counter = zeros(Int, 9)
    
    # Loop through all cells and calculate packing densities, then concentration of species 1
    for i in 1:num_particles
        x_idx, y_idx, z_idx = _compute_cell_index_cartesian(x_data[i], y_data[i], z_data[i], x_min, y_min, z_min, recip_dx, recip_dy, recip_dz)
        base = [x_idx, y_idx, z_idx]
        cell_idx = _find_global_cell_index_cartesian(x_idx, y_idx, z_idx, x_divisions, y_divisions)
        overlap_values, overlaps = _particle_overlaps_cartesian(mesh_boundaries[cell_idx, :], x_data[i], y_data[i], z_data[i], radii[i])
        num_overlaps = sum(overlaps)

        overlaps_x = any(overlaps[1:2])
        overlaps_y = any(overlaps[3:4])
        overlaps_z = any(overlaps[5:6])

        if all(overlaps[1:2]) || all(overlaps[3:4]) || all(overlaps[5:6])
            # Double overlap detected, assign full volume to cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        end

        # The below now assumes no two boundaries of the same dimension are overlapped simultaneously

        relative_x, relative_y, relative_z = 0, 0, 0

        if overlaps_x relative_x = overlaps[1] ? -1 : 1 end
        if overlaps_y relative_y = overlaps[3] ? -1 : 1 end
        if overlaps_z relative_z = overlaps[5] ? -1 : 1 end

        relatives = [relative_x, relative_y, relative_z]

        overlap_dims = [i for i in 1:3 if
                                 (i == 1 && overlaps_x) ||
                                 (i == 2 && overlaps_y) ||
                                 (i == 3 && overlaps_z)]

        # case = 0

        if num_overlaps == 0
            # Particle completely inside cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        elseif num_overlaps == 1
            overlap_value = overlap_values[overlaps][1]
            # Overlaps one boundary, giving a single spherical cap
            main_volume = single_cap_intersection(radii[i], overlap_value)
            remainder_volumes = [particle_volumes[i] - main_volume]

            remainder_indices = [[x_idx + relative_x, y_idx + relative_y, z_idx + relative_z]]
        elseif num_overlaps == 2
            overlap_vals = overlap_values[overlaps]
            
            # Compute the main volume for the cell containing the particle's center
            main_volume = double_cap_intersection(radii[i], overlap_vals[1], overlap_vals[2])
            
            # Compute cap volumes for each overlapping boundary individually
            cap_volume_1 = single_cap_intersection(radii[i], overlap_vals[1])
            cap_volume_2 = single_cap_intersection(radii[i], overlap_vals[2])
            
            # Derive remainder volumes ensuring volume conservation
            remainder_volume_1 = cap_volume_2 - main_volume
            remainder_volume_2 = cap_volume_1 - main_volume
            remainder_volume_3 = particle_volumes[i] - (main_volume + remainder_volume_1 + remainder_volume_2)
            
            # Compute the remainder cell indices using relative offsets
            r1 = copy(base)
            r1[overlap_dims[1]] += relatives[overlap_dims[1]]
            
            r2 = copy(base)
            r2[overlap_dims[2]] += relatives[overlap_dims[2]]
            
            r3 = copy(base)
            r3[overlap_dims[1]] += relatives[overlap_dims[1]]
            r3[overlap_dims[2]] += relatives[overlap_dims[2]]
            
            remainder_indices = [r1, r2, r3]
            remainder_volumes = [remainder_volume_1, remainder_volume_2, remainder_volume_3]
        elseif num_overlaps == 3
            # Instead of computing triple-cap intersections, split the particle volume equally.
            main_volume = particle_volumes[i] / 8
        
            # Generate offsets for the 7 remainder cells (exclude the base cell [0,0,0])
            offsets = [
                [relatives[1], 0, 0],
                [0, relatives[2], 0],
                [0, 0, relatives[3]],
                [relatives[1], relatives[2], 0],
                [relatives[1], 0, relatives[3]],
                [0, relatives[2], relatives[3]],
                [relatives[1], relatives[2], relatives[3]]
            ]
        
            # Compute remainder indices by adding each offset to the base index
            remainder_indices = [base .+ off for off in offsets]
        
            # Create a vector of 7 elements, each equal to subvol (1/8 of the particle volume)
            remainder_volumes = fill(main_volume, 7)
        else
            # println("Calculated an invalid number of overlaps: $num_overlaps")
            # println("This is likely due to a cell somewhere being thinner than a particle")
            # Assign particle completely inside cell
            main_volume = particle_volumes[i]
            particle_volume_per_cell[cell_idx] += main_volume
            continue
        end

        # problem_remainder_volumes = findall(remainder_volumes .< -1e-11)
        
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
            if any(remainder_idx .== 0) || any(remainder_idx .> [x_divisions, y_divisions, z_divisions])
                # Remainder cell outside boundaries
                continue
            end
            # println(remainder_idx)
            remainder_global_idx = _find_global_cell_index_cartesian(remainder_idx..., x_divisions, y_divisions)
            particle_volume_per_cell[remainder_global_idx] += remainder_volume
        end

    end
    # println(case_counter)
    # println(total_case_counter)
    return particle_volume_per_cell
end


@inline function _particle_overlaps_cartesian(cell_boundaries, x, y, z, radius)
    x_min, x_max, y_min, y_max, z_min, z_max = cell_boundaries
    tolerance = 1e-10
    overlap_values = [
        x_min - x,
        x - x_max,
        y_min - y,
        y - y_max,
        z_min - z,
        z - z_max
    ]
    
    overlaps = -radius - tolerance .< overlap_values .< radius + tolerance

    
    return overlap_values, overlaps
end


@inline function _find_global_cell_index_cartesian(x_idx, y_idx, z_idx, x_divisions, y_divisions)
    
    global_cell_index = x_idx + (y_idx - 1) * x_divisions + (z_idx - 1) * x_divisions * y_divisions

    return global_cell_index
end


@inline function _compute_cell_index_cartesian(x, y, z, x_min, y_min, z_min, recip_dx, recip_dy, recip_dz)
    # Compute x_idx, y_idx, z_idx
    x_idx = Int(floor((x - x_min) * recip_dx)) + 1
    y_idx = Int(floor((y - y_min) * recip_dy)) + 1
    z_idx = Int(floor((z - z_min) * recip_dz)) + 1

    return x_idx, y_idx, z_idx
end


@inline function _is_particle_outside(x, y, z, boundaries)
    return (x < boundaries.x_min) || (x > boundaries.x_max) ||
           (y < boundaries.y_min) || (y > boundaries.y_max) ||
           (z < boundaries.z_min) || (z > boundaries.z_max)
end


end # module Cartesian