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

        # Ensure min < max for each dimension (except theta)
        if system == :cartesian
            for i in 1:2:length(required_keys)
                min_key, max_key = required_keys[i], required_keys[i + 1]
                if boundaries_input[min_key] >= boundaries_input[max_key]
                    println("Error-causing boundaries: $boundaries_input")
                    throw(ArgumentError("Invalid boundaries: '$min_key' must be less than '$max_key'."))
                end
            end
        elseif system == :cylindrical
            if boundaries_input[:r_min] >= boundaries_input[:r_max]
                throw(ArgumentError("r_min must be less than r_max. r_min: $(boundaries_input[:r_min]), r_max: $(boundaries_input[:r_max])"))
            end
            if boundaries_input[:z_min] >= boundaries_input[:z_max]
                throw(ArgumentError("z_min must be less than z_max. z_min: $(boundaries_input[:z_min]), z_max: $(boundaries_input[:z_max])"))
            end
        else
            throw(ArgumentError("Invalid system: $system. Choose :cartesian or :cylindrical"))
        end

        # Convert dictionary to array
        return Float64[boundaries_input[key] for key in required_keys]

    elseif isa(boundaries_input, AbstractVector)
        # Handle array input
        if length(boundaries_input) != 6
            throw(ArgumentError("Expected a boundary array of length 6, but got length $(length(boundaries_input))."))
        end

        # Ensure min < max for each dimension (except theta)
        if system == :cartesian
            for i in 1:2:5
                if boundaries_input[i] >= boundaries_input[i + 1]
                    println("Error-causing boundaries: $boundaries_input")
                    throw(ArgumentError("Invalid boundaries: boundaries[$i] ($(boundaries_input[i])) must be less than boundaries[$i+1] ($(boundaries_input[i+1]))."))
                end
            end
        elseif system == :cylindrical
            if boundaries_input[1] >= boundaries_input[2]
                throw(ArgumentError("r_min must be less than r_max. r_min: $(boundaries_input[1]), r_max: $(boundaries_input[2])"))
            end
            if boundaries_input[5] >= boundaries_input[6]
                throw(ArgumentError("z_min must be less than z_max. z_min: $(boundaries_input[5]), z_max: $(boundaries_input[6])"))
            end
        else
            throw(ArgumentError("Invalid system: $system. Choose :cartesian or :cylindrical"))
        end

        return Float64.(boundaries_input)  # Ensure the output is Float64

    else
        throw(ArgumentError("Boundaries must be either a dictionary or an array."))
    end
end


end # module Utils