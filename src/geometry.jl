# geometry.jl
# This module contains all the necessary geometrical functions required by
# the packing density calculator

module Geometry

export convert_to_cylindrical

function simpson(x, y)
    """
    Numerically integrates a function using Simpson's 1/3 rule for evenly
    spaced points.

    If the number of points is odd, the entire range is integrated using
    Simpson's rule. If the number of points is even, Simpson's rule is applied
    to the first `n-1` points, and the trapezoidal rule is applied to the last
    interval.

    Parameters
    ----------
    x : array_like
        1D array of x-coordinates at which the function is sampled. The
        x-coordinates must be evenly spaced within a tolerance of `1e-10`.
    y : array_like
        1D array of y-coordinates (function values) corresponding to the
        x-coordinates.

    Returns
    -------
    float
        The estimated integral of the function.

    Raises
    ------
    ValueError
        If `x` and `y` have different lengths.
        If `x` is not evenly spaced within a tolerance of `1e-10`.

    Example
    --------
    Integrate a cubic function between 1 and 5:
    
    >>> import numpy as np
    >>> def func(x):
    ...     return x**3 - 6*x**2 + 11*x - 6
    >>> x = np.linspace(1, 5, num=1000)
    >>> y = func(x)
    >>> result = simpsons(x, y)
    16.0000  # Approximation of the integral
    """

    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length."))
    end

    # Check that points are evenly spaced
    spacings = diff(x)
    if !isapprox(spacings, spacings[1]*ones(length(x)-1); atol=1e-7)
        println(spacings)
        throw(ArgumentError("x-coordinates are not evenly spaced within
                         tolerance of epsilon = 1e-7."))
    end

    num_points = length(x)

    if num_points % 2 == 0
        @warn "Warning: Even number of points detected. Applying trapezoidal \
rule for the last interval."
        h = x[2] - x[1]
        result = (h / 3) * (
            y[begin] + y[end-1] + 4 * sum(y[2:2:end-1]) + 2 * sum(y[3:2:end-3])
                           )
        # Apply trapezoidal rule for the last interval
        result += 0.5 * (x[end] - x[end-1]) * (y[end] + y[end-1])
        return result
    end

    # Standard Simpson's rule
    h = x[2] - x[1]
    result = (h / 3) * (
        y[begin] + y[end] + 4 * sum(y[2:2:end]) + 2 * sum(y[3:2:end-2])
                        )
    return result
end


function triple_cap_analytical(R_p, a, b, z_min, z_max)

    # Predefined constants (using multiplication rather than division)
    half      = 0.5
    quarter   = 0.25
    one_third = 0.3333333333333333
    one_sixth = 0.16666666666666666
    pi_val    = pi

    # Helper functions to safely compute square roots and inverse sines.
    # Instead of throwing an error on a negative radicand we use max(0,·).
    safe_sqrt(x) = sqrt(max(0, x))
    # For division, if the denominator is near zero we return zero.
    safe_div(x, y) = abs(y) < 1e-12 ? 0.0 : x / y
    # Clamp the argument for asin into [-1, 1] to prevent domain errors.
    safe_asin(x) = asin(clamp(x, -1, 1))

    # Precompute quantities that do not depend on z.
    R_p2   = R_p * R_p
    R_p3   = R_p2 * R_p
    sqrt_Rpa = safe_sqrt(R_p2 - a * a)
    sqrt_Rpb = safe_sqrt(R_p2 - b * b)

    # Define a helper function to compute the antiderivative at a given z.
    # (This corresponds to the integrals computed in triple_cap_antiderivative.)
    function F(z)
        # Use safe_sqrt to ensure nonnegative radicands.
        sqrt_Rpz  = safe_sqrt(R_p2 - z*z)
        sqrt_Rpaz = safe_sqrt(R_p2 - a*a - z*z)
        sqrt_Rpbz = safe_sqrt(R_p2 - b*b - z*z)

        # The trigonometric functions in the expression depend on z.
        # We use safe_div to avoid division-by-zero and safe_asin to clamp arguments.
        asin_z_a  = safe_asin(safe_div(z, sqrt_Rpa))
        asin_a_z  = safe_asin(safe_div(a, sqrt_Rpz))
        atan_a    = atan(safe_div(a * z, R_p * (sqrt_Rpaz + 1e-12)))
        asin_z_b  = safe_asin(safe_div(z, sqrt_Rpb))
        asin_b_z  = safe_asin(safe_div(b, sqrt_Rpz))
        atan_b    = atan(safe_div(b * z, R_p * (sqrt_Rpbz + 1e-12)))
        z3        = z * z * z

        T1 = quarter * pi_val * (R_p2 * z - one_third * z3)
        T2 = ((one_sixth * z3 - half * R_p2 * z) * asin_a_z + one_third * R_p3 * atan_a)
        T3 = (one_sixth * a * (a*a - 3 * R_p2) * asin_z_a - one_third * a * z * sqrt_Rpaz)
        T4 = ((one_sixth * z3 - half * R_p2 * z) * asin_b_z + one_third * R_p3 * atan_b)
        T5 = (one_sixth * b * (b*b - 3 * R_p2) * safe_asin(safe_div(b, sqrt_Rpz)) - one_third * b * z * sqrt_Rpbz)
        T6 = a * b * z

        return T1 + T2 + T3 + T4 + T5 + T6
    end

    # Return the definite volume as the difference between the antiderivative evaluated at z_max and z_min.
    return F(z_max) - F(z_min)
end



function triple_cap_integrator(R, a, b,
                          c_lim_lower, c_lim_upper;
                          num_simpson_sample_points=7)
    """
    Function for integrating the differential volume of slices of a double
    spherical cap intersection. R is the radius of the sphere, a and b are the
    distances of two planes from the sphere's centre. c_lim_lower and
    c_lim_upper are the integration limits in the third dimension.
    """

    # 6 sample points, can be improved if you're being very precise
    if isnan(c_lim_lower)
        throw(ArgumentError("c lower contains NaN values."))
    elseif isnan(c_lim_upper)
        throw(ArgumentError("c upper contains NaN values."))
    end
    c_values = LinRange(c_lim_lower, c_lim_upper,
                        num_simpson_sample_points)
    radius_values = sqrt.(max.(0, R^2 .- c_values.^2))
    # cross_sectional_area_values = [double_circle_intersection(r, a, b) for r in radius_values]
    cross_sectional_area_values = double_circle_intersection.(radius_values, a, b)
    # Integrate cross sectional slice throughout the volume
    volume = simpson(c_values, cross_sectional_area_values)
    return volume
end


function double_circle_intersection(r, a, b)
    """
    Calculate the area of intersection between a circle and two perpendicular
    chords.

    Args:
        r (float): Radius of the circle.
        a (float): Distance of the first chord's center from the circle center.
        b (float): Distance of the second chord's center from the circle center

    Returns:
        float: The area of the intersected region.

    Raises:
        ValueError: If the geometry does not permit an intersection.

    :Example:
    >>> double_circle_intersection(2, 0, 0)
    3.141592653589793
    """
    # Evaluating cases
    if a^2 + b^2 < r^2
        # a and b are contained within the circle, double intersection exists

        # Calculate the sector area (the area of the circle wedge)
        sector_area = 0.5 * r^2 * (3*pi/2 + asin(a / r) + asin(b / r))

        # Calculate areas of the two remaining triangles
        triangle_area = 0.5 * a * (sqrt(r^2 - a^2) - b) + 0.5 * b * (sqrt(r^2 - b^2) - a)

        # Circle area minus sector area and triangles
        intersection_area = pi*r^2 - sector_area - triangle_area

        return intersection_area
    else
        # a and b outside or on circle
        # Special cases, a double intersection does not exist
        if (a >= 0 && b >= 0) || a >= r || b >= r
            return 0  # No circle
        elseif a <= -r && b <= -r
            return pi*r^2  # Full circle
        # Segment with centre-chord distance = b
        elseif (a < 0 && 0 < b < r) || (a <= -r && -r < b <= 0)
            return r^2 * (acos(b/r) - (b/r) * sqrt(1 - (b/r)^2))
        # Segment with centre-chord distance = a
        elseif (b < 0 && 0 < a < r) || (b <= -r && -r < a <= 0)
            return r^2 * (acos(a/r) - (a/r) * sqrt(1 - (a/r)^2))
        elseif -r < a < 0 && -r < b < 0
            # Circle missing two minor segments, with centre-chord distances
            # -a and -b
            return pi*r^2 - (
                r^2 * (acos(-b/r) + (b/r) * sqrt(1 - (b/r)^2))) - (
                r^2 * (acos(-a/r) + (a/r) * sqrt(1 - (a/r)^2)))
        end
    end
end


@inline function single_cap_intersection(R, a)
    """
    Function that evaluates the (analytical) volume of a spherical cap. The
    sphere has radius R, and the distance from the sphere's centre to the
    boundary is a.

    :Example:
    >>> single_cap_intersection(1, 0.2)
    1.47445415208481
    """
    return pi * (R - a)^2 * (2 * R + a) / 3
end

@inline function double_cap_intersection(R, a, b)
    """
    Function that evaluates the volume of a double spherical cap intersection.
    The sphere has radius R, and the distances from the centre to the
    boundaries are a and b

    The cross sectional area function can safely just be integrated from -R
    to R. However, this may be wasteful as this can include regions where the
    cross sectional area is zero. The integration limits are set as small as
    possible, such that they just encapsulate the cap volume.

    :Example:
    >>> double_cap_intersection(1, 0.2, 0.3)
    0.3974826065772735
    """
    
    if a^2 + b^2 <= R^2
        # a and b are contained within sphere, double cap intersection exists
        if a < 0 && b < 0
            c_lim_upper = R
        elseif a < 0
            c_lim_upper = sqrt(max(0, R^2 - b^2))
        elseif b < 0
            c_lim_upper = sqrt(max(0, R^2 - a^2))
        else
            c_lim_upper = sqrt(max(0, R^2 - a^2 - b^2))
        end
    else
        # Short-circuiting for cases which have analytical solutions
        # (perfect accuracy and reduces computational load)
        if a > 0 && b > 0
            # No intersection
            return 0
        elseif a < 0 && b > 0
            # Single cap intersection, with centre-chord distance = b
            return pi * (R - b)^2 * (3 * R - (R - b)) / 3
        elseif b < 0 && a > 0
            # Single cap intersection, with centre-chord distance = a
            return pi * (R - a)^2 * (3 * R - (R - a)) / 3
        else
            # Sphere missing two caps, with centre-chord distances -a and -b
            return 4/3 * pi * R^3 - (
                pi * (R + a)^2 * (3 * R - (R + a)) / 3) - (
                pi * (R + b)^2 * (3 * R - (R + b)) / 3)
        end

    end

    # The double cap intersection is symmetrical, so c_lim_lower is set to 0
    # and the volume doubled
    c_lim_lower = 0
    # return 2*triple_cap_integrator(R, a, b, c_lim_lower, c_lim_upper;
    #                                num_simpson_sample_points=3)
    return 2*triple_cap_analytical(R, a, b, c_lim_lower, c_lim_upper)
end


@inline function triple_cap_intersection(R, a, b, c)
    """
    Function that evaluates the volume of a triple cap intersection. The sphere
    has radius R, and the distance from the sphere's centre to the boundaries
    are a, b and c.

    The cross sectional area function must now be carefully integrated to
    include the intersection with the boundary defined by c. The upper
    integration limit is set as low as possible, such that it still entirely
    encapsulates the cap volume. The lower integration limit is set as c,
    unless the cap is symmetrical (c <= -c_lim_upper) or there is no
    intersection (c >= c_lim_upper).

    :Example:
    >>> triple_cap_intersection(1, 0.3, 0.1, 0.2)
    0.16451538109365088
    """

    if a^2 + b^2 <= R^2
        # a and b are contained within sphere
        # This means a triple cap intersection can exist (depending on c)
        if a < 0 && b < 0
            c_lim_upper = R
        elseif a < 0
            c_lim_upper = sqrt(max(0, R^2 - b^2))
        elseif b < 0
            c_lim_upper = sqrt(max(0, R^2 - a^2))
        else
            c_lim_upper = sqrt(max(0, R^2 - a^2 - b^2))
        end
    else
        # Short-circuiting for cases which have analytical solutions
        # (perfect accuracy and reduces computational load)
        if a > 0 && b > 0
            # No intersection
            return 0
        elseif a < 0 && b > 0
            if c <= -sqrt(R^2 - b^2)
                # Single cap intersection, with centre-chord distance = b
                return pi * (R - b)^2 * (3 * R - (R - b)) / 3
            elseif c >= sqrt(R^2 - b^2)
                # No intersection
                return 0
            else
                c_lim_upper = sqrt(max(0, R^2 - b^2))
            end
        elseif b < 0 && a > 0
            if c <= -sqrt(max(0, R^2 - a^2))
                # Single cap intersection, with centre-chord distance = a
                return pi * (R - a)^2 * (3 * R - (R - a)) / 3
            elseif c >= sqrt(max(0, R^2 - a^2))
                # No intersection
                return 0
            else
                c_lim_upper = sqrt(max(0, R^2 - a^2))
            end
        elseif c > 0 && a < -sqrt(max(0, R^2 - c^2)) && b < -sqrt(max(0, R^2 - c^2))
            # Single cap intersection, with centre-chord distance = c
            return pi * (R - c)^2 * (3 * R - (R - c)) / 3
        elseif b < 0 && a < 0
            if c <= -max(sqrt(max(0, R^2 - a^2)), sqrt(max(0, R^2 - b^2)))
                # Sphere missing three single caps, with centre-chord distances
                # -a, -b, and -c
                return 4/3 * pi * R^3 - (
                    pi * (R + a)^2 * (3 * R - (R + a)) / 3) - (
                    pi * (R + b)^2 * (3 * R - (R + b)) / 3) - (
                    pi * (R + c)^2 * (3 * R - (R + c)) / 3)
            else
                c_lim_upper = R
            end
        else
            c_lim_upper = R
        end
    end

    if c >= c_lim_upper
        # No intersection
        return 0
    elseif c <= -c_lim_upper
        # Symmetrical -> double cap intersection
        c_lim_lower = -c_lim_upper
    else
        # c intersects the double cap intersection
        # -> integrate between c and c_lim_upper
        c_lim_lower = c
    end

    # return triple_cap_integrator(R, a, b, c_lim_lower, c_lim_upper)
    return triple_cap_analytical(R, a, b, c_lim_lower, c_lim_upper)

end


function circle_rectangle_intersection(x_min, x_max, y_min, y_max, circle_radius)
    """
    Compute the area of intersection between a circle and a rectangle.

    Args:
        x_min::Float64: Minimum x-coordinate of the rectangle.
        x_max::Float64: Maximum x-coordinate of the rectangle.
        y_min::Float64: Minimum y-coordinate of the rectangle.
        y_max::Float64: Maximum y-coordinate of the rectangle.
        circle_radius::Float64: Radius of the circle, assumed to be centered at
                                the origin (0, 0).

    Returns:
        Float64: The area of the intersection between the circle and the
                 rectangle.
    """
    # Compute distances of all corners from the circle's center (x=0, y=0)
    corners = [
        (x_min, y_min), (x_min, y_max),
        (x_max, y_min), (x_max, y_max)
    ]
    distances = [sqrt(x^2 + y^2) for (x, y) in corners]

    num_outside = sum(d >= circle_radius for d in distances)

    # Determine closest and furthest x and y coordinates
    if abs(x_min) <= abs(x_max)
        x_closest, x_furthest = x_min, x_max
        x_inv = false
    else
        x_closest, x_furthest = x_max, x_min
        x_inv = true
    end

    if abs(y_min) <= abs(y_max)
        y_closest, y_furthest = y_min, y_max
        y_inv = false
    else
        y_closest, y_furthest = y_max, y_min
        y_inv = true
    end

    # Match-case equivalent in Julia using a conditional block
    if num_outside == 0
        # Fully inside the circle
        return (x_max - x_min) * (y_max - y_min)

    elseif num_outside == 4
        # Fully outside the circle
        critical_boundary = sort([abs(x_min), abs(x_max), abs(y_min), abs(y_max)])[3]

        straddling_axis = (x_min * x_max < 0) || (y_min * y_max < 0)

        if critical_boundary < circle_radius && straddling_axis
            overlap_factor = critical_boundary / circle_radius
            overlap_area = circle_radius^2 * (acos(overlap_factor) -
                                              (overlap_factor) * sqrt(1 - (overlap_factor)^2))
            return overlap_area
        else
            return 0
        end

    elseif num_outside == 3
        # Exactly three corners are outside
        a = ifelse(x_inv, -x_closest, x_closest)
        b = ifelse(y_inv, -y_closest, y_closest)
        return double_circle_intersection(circle_radius, a, b)

    elseif num_outside == 2
        # Exactly two corners are outside
        a = ifelse(x_inv, -x_closest, x_closest)
        b = ifelse(y_inv, -y_closest, y_closest)
        intersection_area = double_circle_intersection(circle_radius, a, b)

        # Additional correction for partial overlap
        a_extra1 = max(ifelse(x_inv, -x_closest, x_closest),
                       ifelse(y_inv, -y_closest, y_closest))
        b_extra1 = min(ifelse(y_inv, -y_furthest, y_furthest),
                       ifelse(x_inv, -x_furthest, x_furthest))

        extra_intersection_area = double_circle_intersection(circle_radius, a_extra1, b_extra1)

        furthest_boundary = max(abs(x_min), abs(x_max), abs(y_min), abs(y_max))
        straddling_axis = (x_min * x_max < 0) || (y_min * y_max < 0)

        overlap_area = if furthest_boundary < circle_radius && straddling_axis
            overlap_factor = furthest_boundary / circle_radius
            circle_radius^2 * (acos(overlap_factor) - overlap_factor * sqrt(1 - overlap_factor^2))
        else
            0
        end

        return intersection_area - extra_intersection_area - overlap_area

    elseif num_outside == 1
        # Exactly one corner is outside
        a = ifelse(x_inv, -x_closest, x_closest)
        b = ifelse(y_inv, -y_closest, y_closest)
        intersection_area = double_circle_intersection(circle_radius, a, b)

        # Additional correction for partial overlap
        a_extra1 = max(ifelse(x_inv, -x_closest, x_closest),
                       ifelse(y_inv, -y_closest, y_closest))
        b_extra1 = min(ifelse(y_inv, -y_furthest, y_furthest),
                       ifelse(x_inv, -x_furthest, x_furthest))

        a_extra2 = min(ifelse(x_inv, -x_closest, x_closest),
                       ifelse(y_inv, -y_closest, y_closest))
        b_extra2 = max(ifelse(y_inv, -y_furthest, y_furthest),
                       ifelse(x_inv, -x_furthest, x_furthest))

        extra_intersection_area = double_circle_intersection(circle_radius, a_extra1, b_extra1) +
                                  double_circle_intersection(circle_radius, a_extra2, b_extra2)

        return intersection_area - extra_intersection_area

    else
        throw(ArgumentError("Unexpected number of corners outside the circle: $num_outside"))
    end
end


function circle_circle_intersection(R_C, R_c, r)
    """
    Compute the intersection area between two circles in 2D.

    Args:
        R_C::Float64: Radius of the first circle.
        R_c::Float64: Radius of the second circle.
        r::Float64: Distance between the centers of the circles.

    Returns:
        Float64: The area of intersection between the two circles.
    """
    if r == 0
        # Circles are concentric; return the area of the smaller circle
        return pi * min(R_C^2, R_c^2)
    end

    # Compute the distances of the intersection points on the line joining the centres
    r_i = r / 2 + 0.5 * (R_C^2 - R_c^2) / r
    r_ic = r - r_i

    # Use clamp for FP precision errors
    ratio_I1 = clamp(r_ic / R_c, -1, 1)
    ratio_I2 = clamp(r_i / R_C, -1, 1)

    # Calculate the intersection areas for the two circular segments
    I1 = R_c^2 * (acos(ratio_I1) - (ratio_I1) * sqrt(1 - (ratio_I1)^2))
    I2 = R_C^2 * (acos(ratio_I2) - (ratio_I2) * sqrt(1 - (ratio_I2)^2))

    return I1 + I2
end


function circle_annulus_intersection(r_p, r, r_overlap; min_boundary=false)
    """
    Compute the intersection area between a circle and an annular region.

    Args:
        r_p::Float64: Radius of the circle.
        r::Float64: Radius of the annulus.
        r_overlap::Float64: Overlap distance in the radial direction.
        min_boundary::Bool: Whether the overlap is with the minimum radius
                            boundary (default: false).

    Returns:
        Float64: The intersection area between the circle and the annulus.
    """

    # Handle special cases
    if r_overlap >= r_p
        # No intersection
        return 0.0
    elseif r_overlap <= -r_p
        # Full intersection
        return pi * r_p^2
    elseif min_boundary
        # Intersection with the inner boundary of the annulus
        return pi * r_p^2 - circle_circle_intersection(r + r_overlap, r_p, r)
    else
        # Intersection with the outer boundary of the annulus
        return circle_circle_intersection(r - r_overlap, r_p, r)
    end
end


function sphere_cylinder_integrator(R_p, r, r_overlap, z_lim_lower, z_lim_upper;
                                    min_boundary=false, num_simpson_sample_points=7)
    """
    Numerically integrate the volume of a sphere intersecting with a
    cylindrical boundary.

    Args:
        R_p::Float64: Radius of the sphere.
        r::Float64: Radius of the cylinder.
        r_overlap::Float64: Overlap distance in the radial direction.
        z_lim_lower::Float64: Lower limit of integration in the z-direction.
        z_lim_upper::Float64: Upper limit of integration in the z-direction.
        min_boundary::Bool: Whether the overlap is with the minimum radial
                            boundary (default: false).
        num_simpson_sample_points::Int: Number of points for Simpson's rule
                                        integration (default: 7).

    Returns:
        Float64: The volume of the sphere-cylinder intersection.
    """
    if isnan(z_lim_lower)
        throw(ArgumentError("z lower contains NaN values."))
    elseif isnan(z_lim_upper)
        throw(ArgumentError("z upper contains NaN values."))
    end
    # Generate evenly spaced points in the z-direction
    z_values = LinRange(z_lim_lower, z_lim_upper, num_simpson_sample_points)
    
    # Compute the radius of the sphere's cross-sections at each z-value
    r_p_values = sqrt.(max.(0, R_p^2 .- z_values.^2))
    
    # Compute the cross-sectional areas for each slice
    cross_sectional_area_values = [
        circle_annulus_intersection(r_p, r, r_overlap; min_boundary=min_boundary)
        for r_p in r_p_values
    ]

    # Integrate the cross-sectional areas to compute the volume
    volume = simpson(z_values, cross_sectional_area_values)
    return volume
end


function sphere_cylinder_intersection(R_p, r, r_overlap; min_boundary=false)
    """
    Compute the volume of intersection between a sphere and a cylindrical
    boundary.

    Args:
        R_p::Float64: Radius of the sphere.
        r::Float64: Radius of the cylinder.
        r_overlap::Float64: Overlap in the radial direction.
        min_boundary::Bool: Whether the overlap is with the minimum boundary
                            (default: false).

    Returns:
        Float64: Volume of the intersection between the sphere and the cylinder.
    """

    # Determine the integration limits in the z-direction
    if r_overlap >= R_p
        # No intersection
        return 0.0
    elseif r_overlap <= -R_p
        # Full sphere
        return (4 / 3) * pi * R_p^3
    elseif r_overlap > 0
        z_lim_upper = sqrt(max(0, R_p^2 - r_overlap^2))
    else
        z_lim_upper = R_p
    end

    # Symmetrical about z = 0
    z_lim_lower = 0.0

    # Integrate the volume using sphere_cylinder_integrator
    volume = 2 * sphere_cylinder_integrator(R_p, r, r_overlap, z_lim_lower, z_lim_upper;
                                            min_boundary=min_boundary, num_simpson_sample_points=7)
    return volume
end


function sphere_cylinder_plane_intersection(R_p, r, r_overlap, z_overlap; min_boundary=false)
    """
    Calculate the volume of intersection between a sphere, a cylinder, and a plane.

    Args:
        R_p::Float64: Radius of the sphere.
        r::Float64: Radius of the cylinder.
        r_overlap::Float64: Overlap in the radial direction.
        z_overlap::Float64: Overlap in the vertical (z) direction.
        min_boundary::Bool: Whether the overlap is with the minimum boundary
                            (default: false).

    Returns:
        Float64: Volume of the intersection between the sphere, cylinder, and plane.
    """

    # Determine z limits based on r_overlap
    if r_overlap >= R_p
        # No intersection
        return 0.0
    elseif r_overlap <= -R_p
        # Entire sphere intersects the cylinder, but constrained by the plane
        # This is a spherical cap
        return single_cap_intersection(R_p, z_overlap)
    elseif r_overlap > 0
        z_lim_upper = sqrt(max(0, R_p^2 - r_overlap^2))
    else
        z_lim_upper = R_p
    end

    # Check if z_overlap prevents intersection
    if z_overlap > z_lim_upper
        # Plane is above the intersection region; no intersection
        return 0.0
    elseif z_overlap < -z_lim_upper
        # Plane intersects the entire sphere-cylinder volume symmetrically
        z_lim_lower = -z_lim_upper
    else
        # Plane intersects part of the sphere-cylinder intersection
        if isnan(z_overlap)
            print("z overlap NaN")
        end
        z_lim_lower = z_overlap
    end

    # Integrate the volume from z_lim_lower to z_lim_upper
    volume = sphere_cylinder_integrator(R_p, r, r_overlap, z_lim_lower, z_lim_upper;
                                        min_boundary=min_boundary, num_simpson_sample_points=7)
    return volume
end


function compute_cell_volume(; boundaries=nothing, system=:cartesian, cylinder_radius=nothing)
    """
    Calculate the volume of a cuboidal or cylindrical cell, adjusting for
    overlap with a cylindrical boundary when specified for cuboids.

    Args:
        boundaries::Tuple: A tuple defining the cell's boundaries. For :cartesian,
                           (x_min, x_max, y_min, y_max, z_min, z_max). For :cylindrical,
                           (r_min, r_max, theta_min, theta_max, z_min, z_max).
        system::Symbol: Coordinate system type (:cartesian or :cylindrical). Default: :cartesian.
        cylinder_radius::Union{Nothing, Float64}: Radius of the cylinder for overlap
                                                  calculations. Default: nothing.

    Returns:
        Float64: The volume of the cell.
    """

    if system == :cartesian
        # Cartesian cell type
        x_min, x_max, y_min, y_max, z_min, z_max = boundaries

        if isnothing(cylinder_radius)
            # Regular Cartesian volume if no cylinder is defined
            return (x_max - x_min) * (y_max - y_min) * (z_max - z_min)
        else
            # Adjust volume for cylinder intersection in the x-y plane
            base_area = circle_rectangle_intersection(x_min, x_max, y_min, y_max, cylinder_radius)
            return base_area * (z_max - z_min)
        end

    elseif system == :cylindrical
        # Cylindrical cell type
        r_min, r_max, theta_min, theta_max, z_min, z_max = boundaries

        # Compute angular difference
        delta_theta = angular_difference(theta_min, theta_max)

        if r_min < 0
            # Full cylindrical cell
            return pi * r_max^2 * (z_max - z_min)
        elseif theta_min == 0 && theta_max == 2pi
            # Full angular range
            return pi * (r_max^2 - r_min^2) * (z_max - z_min)
        elseif delta_theta > 2pi
            # Handle wraparound cases for angular range > 2pi
            return pi * (r_max^2 - r_min^2) * (z_max - z_min)
        else
            # Partial cylindrical cell
            outer_volume = 0.5 * r_max^2 * delta_theta * (z_max - z_min)
            inner_volume = 0.5 * r_min^2 * delta_theta * (z_max - z_min)
            return outer_volume - inner_volume
        end

    else
        throw(ArgumentError("Invalid coordinate system specified. Use 'cartesian' or 'cylindrical'."))
    end
end


function angular_difference(theta1, theta2)
    """
    Calculate the shortest angular difference between two angles,
    wrapped to the range [0, 2π], treating 0 and 2π as a full-circle difference.

    Args:
        theta1::Union{Float64, AbstractArray}: First angle(s) in radians.
        theta2::Union{Float64, AbstractArray}: Second angle(s) in radians.

    Returns:
        Union{Float64, AbstractArray}: The angular difference(s), wrapped to [0, 2π],
                                       with special handling for 0 and 2π.
    """

    # Compute the difference modulo 2 * pi, broadcasting if needed
    delta_theta = mod.(theta2 .- theta1, 2 * pi)

    # Handle the special case where the difference is 0 but the angles are different
    if isa(theta1, AbstractArray) || isa(theta2, AbstractArray)
        # Arrays: Check element-wise and handle the special case
        full_circle_mask = (theta1 .!= theta2) .& (delta_theta .== 0)
        delta_theta[full_circle_mask] .= 2 * pi
    elseif delta_theta == 0 && theta1 != theta2
        # Scalars: Handle the special case
        delta_theta = 2 * pi
    end

    return delta_theta
end


function theta_within_range(theta_data, theta_min, theta_max, factor; mode="inside")
    """
    Check whether angles are within or outside a periodic range.

    Args:
        theta_data::Union{Float64, AbstractArray}: Array or scalar of angles in radians.
        theta_min::Float64: Minimum bound of the range.
        theta_max::Float64: Maximum bound of the range.
        factor::Float64: Adjustment factor for angular overlap.
        mode::String: "inside" to check within range, "outside" to check outside range.

    Returns:
        Union{Bool, AbstractArray{Bool}}: Boolean(s) indicating whether each angle
                                          is within or outside the range.
    """

    # Adjust bounds for particle radii
    adjusted_min = mod(theta_min + factor, 2 * pi)
    adjusted_max = mod(theta_max - factor, 2 * pi)

    # Create masks for standard and wrapped ranges
    standard_mask = adjusted_min <= adjusted_max
    wrapped_mask = !standard_mask

    if mode == "inside"
        # Inside the range
        if isa(theta_data, AbstractArray)
            standard_range = (theta_data .>= adjusted_min) .& (theta_data .<= adjusted_max)
            wrapped_range = (theta_data .>= adjusted_min) .| (theta_data .<= adjusted_max)
        else
            standard_range = (theta_data >= adjusted_min) && (theta_data <= adjusted_max)
            wrapped_range = (theta_data >= adjusted_min) || (theta_data <= adjusted_max)
        end
    elseif mode == "outside"
        # Outside the range
        if isa(theta_data, AbstractArray)
            standard_range = (theta_data .< adjusted_min) .| (theta_data .> adjusted_max)
            wrapped_range = (theta_data .< adjusted_min) .& (theta_data .> adjusted_max)
        else
            standard_range = (theta_data < adjusted_min) || (theta_data > adjusted_max)
            wrapped_range = (theta_data < adjusted_min) && (theta_data > adjusted_max)
        end
    else
        throw(ArgumentError("Mode must be 'inside' or 'outside'"))
    end

    # Combine results based on range type
    if isa(theta_data, AbstractArray)
        return (standard_mask && standard_range) .| (wrapped_mask && wrapped_range)
    else
        return (standard_mask && standard_range) || (wrapped_mask && wrapped_range)
    end
end


function calculate_angular_overlap_factor(r_data, radii)
    """
    Calculate the angular factor for cylindrical overlap calculations.

    Args:
        r_data::Union{Float64, AbstractArray}: Radial distance(s) of particles from the origin.
        radii::Union{Float64, AbstractArray}: Radius/radii of the particles.

    Returns:
        Union{Float64, AbstractArray}: Angular factor(s) for each particle,
                                       based on the ratio of particle radius to radial distance.
    """

    # Ensure radial distances are safe (prevent division by zero)
    safe_r_data = max.(r_data, radii)

    # Calculate angular factor using arcsin, ensuring the ratio is clipped to [-1, 1]
    factor = asin.(clamp.(radii ./ safe_r_data, -1.0, 1.0))

    return factor
end


function convert_to_cylindrical(x_data, y_data)
    """
    Convert Cartesian coordinates to cylindrical coordinates with theta in [0, 2π].

    Args:
        x_data::Union{Float64, AbstractArray}: x-coordinate(s).
        y_data::Union{Float64, AbstractArray}: y-coordinate(s).

    Returns:
        Tuple{Union{Float64, AbstractArray}, Union{Float64, AbstractArray}}:
            - r_data: Radial distance(s) from the origin.
            - theta_data: Angle(s) in radians from the x-axis, in the range [0, 2π].
    """

    # Compute radial distance
    r_data = sqrt.(x_data .^ 2 .+ y_data .^ 2)

    # Compute angle using atan2, ensuring the result is in [-π, π]
    theta_data = atan.(y_data, x_data)

    # Map theta to the range [0, 2π]
    theta_data = mod.(theta_data .+ 2 * pi, 2 * pi)

    return r_data, theta_data
end

end # module Geometry