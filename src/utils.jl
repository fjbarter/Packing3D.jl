# utils.jl
# This module contains general utilities required by functions in the
# calculator

# is_inside_boundaries and is_outside_boundaries are sufficiently fast, as
# the string comparison is only performed once, and the inside and outside
# masks are evaluated for the entire particle dataset in one go.

module Utils

# include("geometry.jl")

using ..Geometry: angular_difference

export convert_boundaries_dictionary

"""
    compute_automatic_boundaries(x_data::AbstractVector{<:Real}, 
                                 y_data::AbstractVector{<:Real}, 
                                 z_data::AbstractVector{<:Real}, 
                                 r_data::AbstractVector{<:Real}, 
                                 theta_data::AbstractVector{<:Real}, 
                                 system::Symbol; 
                                 padding_factor::Float64=0.1)

Compute the minimum and maximum boundaries for Cartesian or cylindrical coordinate systems.

# Arguments
- `x_data::AbstractVector{<:Real}`: x-coordinates of particles (Cartesian).
- `y_data::AbstractVector{<:Real}`: y-coordinates of particles (Cartesian).
- `z_data::AbstractVector{<:Real}`: z-coordinates of particles.
- `r_data::AbstractVector{<:Real}`: Radial positions of particles (Cylindrical).
- `theta_data::AbstractVector{<:Real}`: Angular positions of particles (Cylindrical).
- `system::Symbol`: `:cartesian` or `:cylindrical`.
- `padding_factor::Float64`: Fraction to shrink or expand boundaries. Defaults to `0.1`.

# Returns
- `Dict{Symbol, Any}`: Dictionary containing the boundaries.

# Raises
- `ArgumentError`: If required data for the specified system is missing.
"""
function compute_automatic_boundaries(x_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                       y_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                       z_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                       r_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                       theta_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                       system::Symbol=:cartesian; 
                                       padding_factor::Float64=0.1)

    if system == :cartesian
        if isnothing(x_data) || isnothing(y_data) || isnothing(z_data)
            throw(ArgumentError("x_data, y_data, and z_data are required for Cartesian boundaries."))
        end

        x_range = maximum(x_data) - minimum(x_data)
        y_range = maximum(y_data) - minimum(y_data)
        z_range = maximum(z_data) - minimum(z_data)

        boundaries = Dict(
            :x_min => minimum(x_data) + padding_factor * x_range,
            :x_max => maximum(x_data) - padding_factor * x_range,
            :y_min => minimum(y_data) + padding_factor * y_range,
            :y_max => maximum(y_data) - padding_factor * y_range,
            :z_min => minimum(z_data) + padding_factor * z_range,
            :z_max => maximum(z_data) - padding_factor * z_range
        )
        return boundaries

    elseif system == :cylindrical
        if isnothing(r_data) || isnothing(z_data)
            throw(ArgumentError("r_data and z_data are required for Cylindrical boundaries."))
        end

        z_range = maximum(z_data) - minimum(z_data)

        boundaries = Dict(
            :r_min => maximum(r_data) * (1 - padding_factor),
            :r_max => -maximum(r_data) * (1 - padding_factor),
            :theta_min => -π,
            :theta_max => 3π,
            :z_min => minimum(z_data) + padding_factor * z_range,
            :z_max => maximum(z_data) - padding_factor * z_range
        )
        return boundaries

    else
        throw(ArgumentError("Invalid system specified. Choose 'cartesian' or 'cylindrical'."))
    end
end


"""
    is_inside_boundaries(x_data::Union{AbstractVector{<:Real}, Nothing}, 
                         y_data::Union{AbstractVector{<:Real}, Nothing}, 
                         z_data::Union{AbstractVector{<:Real}, Nothing},
                         r_data::Union{AbstractVector{<:Real}, Nothing}, 
                         theta_data::Union{AbstractVector{<:Real}, Nothing}, 
                         boundaries::AbstractVector{<:Real}, 
                         radii::AbstractVector{<:Real},
                         factor::Union{Real, Nothing}, 
                         system::Symbol)

Determine which particles are completely inside the defined boundaries
for either Cartesian or cylindrical systems.

# Arguments
- `x_data::Union{AbstractVector{<:Real}, Nothing}`: x-coordinates (Cartesian).
- `y_data::Union{AbstractVector{<:Real}, Nothing}`: y-coordinates (Cartesian).
- `z_data::Union{AbstractVector{<:Real}, Nothing}`: z-coordinates (Cartesian).
- `r_data::Union{AbstractVector{<:Real}, Nothing}`: Radial distances (Cylindrical).
- `theta_data::Union{AbstractVector{<:Real}, Nothing}`: Angular positions (Cylindrical).
- `boundaries::AbstractVector{<:Real}`: Boundary values (e.g., `[x_min, x_max, y_min, y_max, z_min, z_max]`).
- `radii::AbstractVector{<:Real}`: Radii of the particles.
- `factor::Union{Real, Nothing}`: Adjustment factor for angular overlaps (Cylindrical).
- `system::Symbol`: `:cartesian` or `:cylindrical`.

# Returns
- `Vector{Bool}`: Boolean mask indicating whether each particle is inside the boundaries.

# Raises
- `ArgumentError`: If required inputs for the chosen system are missing.
"""
function is_inside_boundaries(; x_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                y_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                z_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                r_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                theta_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                boundaries::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                radii::AbstractVector{<:Real}=nothing,
                                factor::Union{AbstractVector{<:Real}, <:Real, Nothing}=nothing,
                                system::Symbol=:cartesian)

    if system == :cartesian
        if isnothing(x_data) || isnothing(y_data) || isnothing(z_data)
            throw(ArgumentError("x_data, y_data, and z_data are required for Cartesian boundaries."))
        end

        # Extract boundaries
        x_min, x_max, y_min, y_max, z_min, z_max = boundaries

        # Compute boolean mask
        return (
            (x_data .>= x_min .+ radii) .&
            (x_data .<= x_max .- radii) .&
            (y_data .>= y_min .+ radii) .&
            (y_data .<= y_max .- radii) .&
            (z_data .>= z_min .+ radii) .&
            (z_data .<= z_max .- radii)
        )

    elseif system == :cylindrical
        if isnothing(r_data) || isnothing(z_data) || isnothing(theta_data)
            throw(ArgumentError("r_data, theta_data, and z_data are required for Cylindrical boundaries."))
        end

        # Extract boundaries
        r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

        # Full cylinder: no r_min constraint
        if r_min < 0
            return (
                (r_data .<= r_max .- radii) .&
                (z_data .>= z_min .+ radii) .&
                (z_data .<= z_max .- radii)
            )
        end

        # Full ring: no angular range constraints
        if theta_min == 0 && theta_max == 2π
            return (
                (r_data .>= r_min .+ radii) .&
                (r_data .<= r_max .- radii) .&
                (z_data .>= z_min .+ radii) .&
                (z_data .<= z_max .- radii)
            )
        end

        # Handle angular constraints
        if isnothing(factor)
            throw(ArgumentError("factor is required for Cylindrical boundaries with angular constraints."))
        end

        theta_min = mod.(theta_min .+ factor, 2π)
        theta_max = mod.(theta_max .- factor, 2π)

        # Theta condition for periodic ranges
        standard_range = (theta_min <= theta_max)
        theta_inside = ifelse(
            standard_range,
            (theta_data .>= theta_min) .& (theta_data .<= theta_max),
            (theta_data .>= theta_min) .| (theta_data .<= theta_max)
        )

        return (
            (r_data .>= r_min .+ radii) .&
            (r_data .<= r_max .- radii) .&
            theta_inside .&
            (z_data .>= z_min .+ radii) .&
            (z_data .<= z_max .- radii)
        )

    else
        throw(ArgumentError("Invalid system specified. Choose 'cartesian' or 'cylindrical'."))
    end
end


"""
    is_outside_boundaries(x_data::Union{AbstractVector{<:Real}, Nothing}, 
                          y_data::Union{AbstractVector{<:Real}, Nothing}, 
                          z_data::Union{AbstractVector{<:Real}, Nothing},
                          r_data::Union{AbstractVector{<:Real}, Nothing}, 
                          theta_data::Union{AbstractVector{<:Real}, Nothing}, 
                          boundaries::AbstractVector{<:Real}, 
                          radii::AbstractVector{<:Real},
                          factor::Union{Real, Nothing}, 
                          system::Symbol)

Determine which particles are completely outside the defined boundaries
for either Cartesian or cylindrical systems.

# Arguments
- `x_data::Union{AbstractVector{<:Real}, Nothing}`: x-coordinates (Cartesian).
- `y_data::Union{AbstractVector{<:Real}, Nothing}`: y-coordinates (Cartesian).
- `z_data::Union{AbstractVector{<:Real}, Nothing}`: z-coordinates (Cartesian).
- `r_data::Union{AbstractVector{<:Real}, Nothing}`: Radial distances (Cylindrical).
- `theta_data::Union{AbstractVector{<:Real}, Nothing}`: Angular positions (Cylindrical).
- `boundaries::AbstractVector{<:Real}`: Boundary values (e.g., `[x_min, x_max, y_min, y_max, z_min, z_max]`).
- `radii::AbstractVector{<:Real}`: Radii of the particles.
- `factor::Union{Real, Nothing}`: Adjustment factor for angular overlaps (Cylindrical).
- `system::Symbol`: `:cartesian` or `:cylindrical`.

# Returns
- `Vector{Bool}`: Boolean mask indicating whether each particle is outside the boundaries.

# Raises
- `ArgumentError`: If required inputs for the chosen system are missing.
"""
function is_outside_boundaries(; x_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                 y_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                 z_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                 r_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                 theta_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                 boundaries::AbstractVector{<:Real}=nothing,
                                 radii::AbstractVector{<:Real}=nothing,
                                 factor::Union{AbstractVector{<:Real}, <:Real, Nothing}=nothing,
                                 system::Symbol=:cartesian)

    if system == :cartesian
        if isnothing(x_data) || isnothing(y_data) || isnothing(z_data)
            throw(ArgumentError("x_data, y_data, and z_data are required for Cartesian boundaries."))
        end

        # Extract boundaries
        x_min, x_max, y_min, y_max, z_min, z_max = boundaries

        # Compute boolean mask
        return (
            (x_data .<= x_min .- radii) .|
            (x_data .>= x_max .+ radii) .|
            (y_data .<= y_min .- radii) .|
            (y_data .>= y_max .+ radii) .|
            (z_data .<= z_min .- radii) .|
            (z_data .>= z_max .+ radii)
        )

    elseif system == :cylindrical
        if isnothing(r_data) || isnothing(z_data) || isnothing(theta_data)
            throw(ArgumentError("r_data, theta_data, and z_data are required for Cylindrical boundaries."))
        end

        # Extract boundaries
        r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

        # Full cylinder: no r_min constraint
        if r_min < 0
            return (
                (r_data .>= r_max .+ radii) .|
                (z_data .<= z_min .- radii) .|
                (z_data .>= z_max .+ radii)
            )
        end

        # Full ring: no angular range constraints
        if theta_min == 0 && theta_max == 2π
            return (
                (r_data .<= r_min .- radii) .|
                (r_data .>= r_max .+ radii) .|
                (z_data .<= z_min .- radii) .|
                (z_data .>= z_max .+ radii)
            )
        end

        # Handle angular constraints
        if isnothing(factor)
            throw(ArgumentError("factor is required for Cylindrical boundaries with angular constraints."))
        end

        theta_min = mod.(theta_min .- factor, 2π)
        theta_max = mod.(theta_max .+ factor, 2π)

        # Theta condition for periodic ranges
        standard_range = (theta_min <= theta_max)
        theta_outside = ifelse(
            standard_range,
            (theta_data .<= theta_min) .| (theta_data .>= theta_max),
            (theta_data .<= theta_min) .& (theta_data .>= theta_max)
        )

        return (
            (r_data .<= r_min .- radii) .|
            (r_data .>= r_max .+ radii) .|
            theta_outside .|
            (z_data .<= z_min .- radii) .|
            (z_data .>= z_max .+ radii)
        )

    else
        throw(ArgumentError("Invalid system specified. Choose 'cartesian' or 'cylindrical'."))
    end
end


"""
    centre_inside_boundaries(x_data::Union{AbstractVector{<:Real}, Nothing}, 
                             y_data::Union{AbstractVector{<:Real}, Nothing}, 
                             z_data::Union{AbstractVector{<:Real}, Nothing},
                             r_data::Union{AbstractVector{<:Real}, Nothing}, 
                             theta_data::Union{AbstractVector{<:Real}, Nothing},
                             boundaries::AbstractVector{<:Real},
                             system::Symbol)

Determine whether the centre of each particle is inside the defined boundaries
for either Cartesian or Cylindrical coordinate systems.

# Arguments
- `x_data::Union{AbstractVector{<:Real}, Nothing}`: x-coordinates (Cartesian).
- `y_data::Union{AbstractVector{<:Real}, Nothing}`: y-coordinates (Cartesian).
- `z_data::Union{AbstractVector{<:Real}, Nothing}`: z-coordinates (Cartesian and Cylindrical).
- `r_data::Union{AbstractVector{<:Real}, Nothing}`: Radial distances (Cylindrical).
- `theta_data::Union{AbstractVector{<:Real}, Nothing}`: Angular positions (Cylindrical).
- `boundaries::AbstractVector{<:Real}`: Boundary values.
- `system::Symbol`: `:cartesian` or `:cylindrical`.

# Returns
- `Vector{Bool}`: Boolean mask indicating whether each particle's centre is inside the boundaries.

# Raises
- `ArgumentError`: If required inputs for the chosen system are missing.
"""
function centre_inside_boundaries(; x_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                   y_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                   z_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                   r_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                   theta_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                   boundaries::AbstractVector{<:Real},
                                   system::Symbol=:cartesian)
    
    if system == :cartesian
        if isnothing(x_data) || isnothing(y_data) || isnothing(z_data)
            throw(ArgumentError("x_data, y_data, and z_data are required for Cartesian boundaries."))
        end

        # Extract boundaries
        x_min, x_max, y_min, y_max, z_min, z_max = boundaries

        # Compute boolean mask
        return (
            (x_data .>= x_min) .&
            (x_data .<= x_max) .&
            (y_data .>= y_min) .&
            (y_data .<= y_max) .&
            (z_data .>= z_min) .&
            (z_data .<= z_max)
        )

    elseif system == :cylindrical
        if isnothing(r_data) || isnothing(theta_data) || isnothing(z_data)
            throw(ArgumentError("r_data, theta_data, and z_data are required for Cylindrical boundaries."))
        end

        # Extract boundaries
        r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

        # Handle angular constraints
        theta_min = mod(theta_min, 2π)
        theta_max = mod(theta_max, 2π)

        standard_range = theta_min <= theta_max
        theta_inside = ifelse(
            standard_range,
            (theta_data .>= theta_min) .& (theta_data .<= theta_max),
            (theta_data .>= theta_min) .| (theta_data .<= theta_max)
        )

        return (
            (r_data .>= r_min) .&
            (r_data .<= r_max) .&
            theta_inside .&
            (z_data .>= z_min) .&
            (z_data .<= z_max)
        )

    else
        throw(ArgumentError("Invalid system specified. Choose :cartesian or :cylindrical."))
    end
end



"""
    calculate_overlaps(x_data::Union{AbstractVector{<:Real}, Nothing},
                       y_data::Union{AbstractVector{<:Real}, Nothing},
                       z_data::Union{AbstractVector{<:Real}, Nothing},
                       r_data::Union{AbstractVector{<:Real}, Nothing},
                       theta_data::Union{AbstractVector{<:Real}, Nothing},
                       radii::AbstractVector{<:Real},
                       boundaries::AbstractVector{<:Real},
                       factor::Union{Real, Nothing},
                       system::Symbol)

Calculate boolean masks for particles overlapping with each boundary
in either Cartesian or cylindrical coordinate systems.

# Arguments
- `x_data::Union{AbstractVector{<:Real}, Nothing}`: x-coordinates (Cartesian).
- `y_data::Union{AbstractVector{<:Real}, Nothing}`: y-coordinates (Cartesian).
- `z_data::Union{AbstractVector{<:Real}, Nothing}`: z-coordinates (Cartesian).
- `r_data::Union{AbstractVector{<:Real}, Nothing}`: Radial distances (Cylindrical).
- `theta_data::Union{AbstractVector{<:Real}, Nothing}`: Angular positions (Cylindrical).
- `radii::AbstractVector{<:Real}`: Radii of the particles.
- `boundaries::AbstractVector{<:Real}`: Boundary values (e.g., `[x_min, x_max, y_min, y_max, z_min, z_max]`).
- `factor::Union{Real, Nothing}`: Adjustment factor for angular overlaps (Cylindrical).
- `system::Symbol`: `:cartesian` or `:cylindrical`.

# Returns
- `Dict{Symbol, AbstractVector{Bool}}`: Boolean masks for overlaps with each boundary.
"""
function calculate_overlaps(;   x_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                y_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                z_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                r_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                theta_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                radii::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                boundaries::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                factor::Union{AbstractVector{<:Real}, <:Real, Nothing}=nothing,
                                system::Symbol=:cartesian)
    
    if system == :cartesian
        if isnothing(x_data) || isnothing(y_data) || isnothing(z_data)
            throw(ArgumentError("x_data, y_data, and z_data are required for Cartesian overlaps."))
        end

        # Extract boundaries
        x_min, x_max, y_min, y_max, z_min, z_max = boundaries

        # Calculate overlaps
        overlaps = Dict(
            :x_min => (x_data .> x_min .- radii) .& (x_data .< x_min .+ radii),
            :x_max => (x_data .> x_max .- radii) .& (x_data .< x_max .+ radii),
            :y_min => (y_data .> y_min .- radii) .& (y_data .< y_min .+ radii),
            :y_max => (y_data .> y_max .- radii) .& (y_data .< y_max .+ radii),
            :z_min => (z_data .> z_min .- radii) .& (z_data .< z_min .+ radii),
            :z_max => (z_data .> z_max .- radii) .& (z_data .< z_max .+ radii)
        )
        return overlaps

    elseif system == :cylindrical
        if isnothing(r_data) || isnothing(z_data) || isnothing(theta_data)
            throw(ArgumentError("r_data, theta_data, and z_data are required for Cylindrical overlaps."))
        end

        # Extract boundaries
        r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

        # Helper function for angular overlaps
        function theta_overlap(theta, t_min, t_max)
            standard_range = (t_min <= t_max)
            return ifelse(
                standard_range,
                (theta .> t_min) .& (theta .< t_max),
                (theta .> t_min) .| (theta .< t_max)
            )
        end

        # Adjust theta boundaries with factor if provided
        if isnothing(factor)
            factor = 0.0  # Default to no angular adjustment
        end
        theta_min_wrapped = mod.(theta_min .- factor, 2π)
        theta_max_wrapped = mod.(theta_max .+ factor, 2π)

        # Handle full cylinder case
        if r_min < 0
            return Dict(
                :r_min => falses(length(r_data)),
                :r_max => (r_data .> r_max .- radii) .& (r_data .< r_max .+ radii),
                :theta_min => falses(length(theta_data)),
                :theta_max => falses(length(theta_data)),
                :z_min => (z_data .> z_min .- radii) .& (z_data .< z_min .+ radii),
                :z_max => (z_data .> z_max .- radii) .& (z_data .< z_max .+ radii)
            )
        end

        # Calculate overlaps for all boundaries
        overlaps = Dict(
            :r_min => (r_data .> r_min .- radii) .& (r_data .< r_min .+ radii),
            :r_max => (r_data .> r_max .- radii) .& (r_data .< r_max .+ radii),
            :theta_min => theta_overlap(theta_data, theta_min_wrapped, theta_min .+ factor),
            :theta_max => theta_overlap(theta_data, theta_max .- factor, theta_max_wrapped),
            :z_min => (z_data .> z_min .- radii) .& (z_data .< z_min .+ radii),
            :z_max => (z_data .> z_max .- radii) .& (z_data .< z_max .+ radii)
        )
        return overlaps

    else
        throw(ArgumentError("Invalid system specified. Choose 'cartesian' or 'cylindrical'."))
    end
end


"""
    calculate_active_overlap_values(total_particles::Int,
                                    x_data::Union{AbstractVector{<:Real}, Nothing},
                                    y_data::Union{AbstractVector{<:Real}, Nothing},
                                    z_data::Union{AbstractVector{<:Real}, Nothing},
                                    r_data::Union{AbstractVector{<:Real}, Nothing},
                                    theta_data::Union{AbstractVector{<:Real}, Nothing},
                                    radii::AbstractVector{<:Real},
                                    boundaries::AbstractVector{<:Real},
                                    overlaps::Dict{Symbol, AbstractVector{Bool}},
                                    system::Symbol)

Calculate the overlap distances for particles intersecting the boundaries
for either Cartesian or cylindrical systems.

# Arguments
- `total_particles::Int`: Total number of particles.
- `x_data::Union{AbstractVector{<:Real}, Nothing}`: x-coordinates (Cartesian).
- `y_data::Union{AbstractVector{<:Real}, Nothing}`: y-coordinates (Cartesian).
- `z_data::Union{AbstractVector{<:Real}, Nothing}`: z-coordinates (Cartesian).
- `r_data::Union{AbstractVector{<:Real}, Nothing}`: Radial distances (Cylindrical).
- `theta_data::Union{AbstractVector{<:Real}, Nothing}`: Angular positions (Cylindrical).
- `radii::AbstractVector{<:Real}`: Radii of the particles.
- `boundaries::AbstractVector{<:Real}`: Boundary values (e.g., `[x_min, x_max, y_min, y_max, z_min, z_max]`).
- `overlaps::Dict{Symbol, AbstractVector{Bool}}`: Boolean masks for overlaps with each boundary.
- `system::Symbol`: `:cartesian` or `:cylindrical`.

# Returns
- `Matrix{Float64}`: A matrix where rows correspond to particles and columns represent boundary overlap distances.
"""
function calculate_active_overlap_values(total_particles::Int;
                                         x_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                         y_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                         z_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                         r_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                         theta_data::Union{AbstractVector{<:Real}, Nothing}=nothing,
                                         boundaries::AbstractVector{<:Real}=nothing,
                                         overlaps::Dict{Symbol, BitVector}=nothing,
                                         system::Symbol=:cartesian)

    # Initialise the result matrix with NaN
    active_overlap_values = fill(NaN, total_particles, 6)

    if system == :cartesian
        if isnothing(x_data) || isnothing(y_data) || isnothing(z_data)
            throw(ArgumentError("x_data, y_data, and z_data are required for Cartesian overlaps."))
        end

        # Extract boundaries
        x_min, x_max, y_min, y_max, z_min, z_max = boundaries

        # Calculate overlap distances for Cartesian boundaries
        active_overlap_values[:, 1] .= ifelse.(overlaps[:x_min], x_min .- x_data, NaN)
        active_overlap_values[:, 2] .= ifelse.(overlaps[:x_max], x_data .- x_max, NaN)
        active_overlap_values[:, 3] .= ifelse.(overlaps[:y_min], y_min .- y_data, NaN)
        active_overlap_values[:, 4] .= ifelse.(overlaps[:y_max], y_data .- y_max, NaN)
        active_overlap_values[:, 5] .= ifelse.(overlaps[:z_min], z_min .- z_data, NaN)
        active_overlap_values[:, 6] .= ifelse.(overlaps[:z_max], z_data .- z_max, NaN)

    elseif system == :cylindrical
        if isnothing(r_data) || isnothing(theta_data) || isnothing(z_data)
            throw(ArgumentError("r_data, theta_data, and z_data are required for Cylindrical overlaps."))
        end

        # Extract boundaries
        r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

        # Radial overlaps
        active_overlap_values[:, 1] .= ifelse.(overlaps[:r_min], r_min .- r_data, NaN)
        active_overlap_values[:, 2] .= ifelse.(overlaps[:r_max], r_data .- r_max, NaN)

        # Angular overlaps
        theta_min_diff = theta_data .- theta_min
        theta_max_diff = theta_max .- theta_data

        # Convert angular differences to distances
        theta_min_dist = r_data .* sin.(mod.(theta_min_diff, 2π))
        theta_max_dist = r_data .* sin.(mod.(theta_max_diff, 2π))

        active_overlap_values[:, 3] .= ifelse.(overlaps[:theta_min], theta_min_dist, NaN)
        active_overlap_values[:, 4] .= ifelse.(overlaps[:theta_max], theta_max_dist, NaN)

        # Axial overlaps
        active_overlap_values[:, 5] .= ifelse.(overlaps[:z_min], z_min .- z_data, NaN)
        active_overlap_values[:, 6] .= ifelse.(overlaps[:z_max], z_data .- z_max, NaN)

    else
        throw(ArgumentError("Invalid system specified. Choose 'cartesian' or 'cylindrical'."))
    end

    return active_overlap_values
end


"""
    convert_boundaries_dictionary(boundaries::Union{Dict{String, <:Real}, AbstractVector{<:Real}},
                                   system::Symbol)

Convert a user-provided boundaries dictionary or validate a boundary array.

# Arguments
- `boundaries::Union{Dict{String, <:Real}, AbstractVector{<:Real}}`:
    - If `Dict`: Specifies boundaries with keys like `"x_min"`, `"x_max"`, etc.
    - If `AbstractVector`: A flat array of 6 elements specifying boundaries.
- `system::Symbol`: The coordinate system, either `:cartesian` or `:cylindrical`.

# Returns
- `AbstractVector{Float64}`: A validated array of boundaries.

# Raises
- `ArgumentError`: If required keys are missing, values are invalid, or the array format is incorrect.
"""
function convert_boundaries_dictionary(boundaries_input::Union{Dict{Symbol, <:Real}, AbstractVector{<:Real}},
                                       system::Symbol)::AbstractVector{Float64}

    if isa(boundaries_input, Dict)
        # Handle dictionary input
        if system == :cartesian
            required_keys = [:x_min, :x_max, :y_min, :y_max, :z_min, :z_max]
        elseif system == :cylindrical
            required_keys = [:r_min, :r_max, :theta_min, :theta_max, :z_min, :z_max]
        else
            throw(ArgumentError("Unsupported system. Choose 'cartesian' or 'cylindrical'."))
        end

        # Check for missing keys
        missing_keys = filter(k -> !haskey(boundaries_input, k), required_keys)
        if !isempty(missing_keys)
            throw(ArgumentError("Missing required boundary keys: $(missing_keys)."))
        end

        # Validate boundary values
        for key in required_keys
            if !(boundaries_input[key] isa Real)
                throw(ArgumentError("Boundary value for '$key' must be a number. Got: $(boundaries_input[key])"))
            end
        end

        # Ensure min < max for each dimension
        for i in 1:2:length(required_keys)
            min_key, max_key = required_keys[i], required_keys[i + 1]
            if boundaries_input[min_key] >= boundaries_input[max_key]
                println("Error-causing boundaries: $boundaries_input")
                throw(ArgumentError("Invalid boundaries: '$min_key' must be less than '$max_key'."))
            end
        end

        # Convert dictionary to array
        return Float64[boundaries_input[key] for key in required_keys]

    elseif isa(boundaries_input, AbstractVector)
        # Handle array input
        if length(boundaries_input) != 6
            throw(ArgumentError("Expected a boundary array of length 6, but got length $(length(boundaries_input))."))
        end

        # Ensure min < max for each dimension
        for i in 1:2:5
            if boundaries_input[i] >= boundaries_input[i + 1]
                println("Error-causing boundaries: $boundaries_input")
                throw(ArgumentError("Invalid boundaries: boundaries[$i] ($(boundaries_input[i])) must be less than boundaries[$i+1] ($(boundaries_input[i+1]))."))
            end
        end

        return Float64.(boundaries_input)  # Ensure the output is Float64

    else
        throw(ArgumentError("Boundaries must be either a dictionary or an array."))
    end
end


end # module Utils