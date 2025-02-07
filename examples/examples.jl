# using Pkg
# Pkg.activate(".")
# base_directory = pwd()
# Pkg.instantiate()
# Pkg.develop(path=joinpath(base_directory, "packing3d_jl_root"))


using Packing3D
using BenchmarkTools
using Plots

function benchmarked_call(f, args...; num_evals=1, num_samples=1, kwargs...)
    result = f(args...; kwargs...)
    benchmark = @benchmark $f($args...; $kwargs...) evals=num_evals samples=num_samples
    return result, benchmark
end

function example1()

    # Define the path to the VTK file
    base_directory = pwd()
    file = joinpath(base_directory, "post", "particles_initial.vtk")
    
    # Check if the file exists
    if !isfile(file)
        error("File not found: $file")
    end

    # Load and process the dataset
    println("Loading file...\n$file")
    data = read_vtk_file(file)
    x_data, y_data, z_data, radii = retrieve_coordinates(data)
    r_data, theta_data = convert_to_cylindrical(x_data, y_data)


    # println(size(radii))

    cylinder_radius = 0.03
    cylinder_base_level = 0.0
    cylinder_height = 0.08

    
    ## PACKING DENSITY ##

    boundaries = Dict(
        :r_min => -1,
        :r_max => cylinder_radius*0.9,
        :theta_min => 0,
        :theta_max => pi,
        :z_min => 0.0025,
        :z_max => 0.0325
    )

    packing_density = compute_packing_cylindrical(; boundaries=boundaries,
                                                    r_data=r_data,
                                                    theta_data=theta_data,
                                                    z_data=z_data,
                                                    radii=radii,
                                                    accurate_cylindrical=true)

    # Split data into two types
    types = data[:point_data][:type]

    type_1 = 1.0
    type_2 = 3.0

    radius_inner = nothing

    small_radius = 0.0005
    large_radius = 0.001

    tolerance = 1e-6
    data_1_mask = abs.(radii .- small_radius) .< tolerance
    data_2_mask = abs.(radii .- large_radius) .< tolerance

    # data_1_mask = abs.(types .- type_1) .< tolerance
    # data_2_mask = abs.(types .- type_2) .< tolerance

    data_1 = extract_points(data, data_1_mask)
    data_2 = extract_points(data, data_2_mask)

    target_num_cells = 200

    SI, SI_benchmark = benchmarked_call(calculate_segregation_intensity,
                                                    data_1,
                                                    data_2;
                                                    cylinder_radius=nothing,
                                                    cylinder_base_level=nothing,
                                                    cylinder_height=cylinder_height,
                                                    target_num_cells=target_num_cells,
                                                    calculate_partial_volumes=true,
                                                    output_num_cells=true,
                                                    clamp_0_to_1=false,
                                                    verbose=false)

    println("Calculated packing density:       $packing_density")
    println("Calculated segregation intensity: $SI")
    println("Elapsed time:                     $(SI_benchmark)")

end

function example2()
    # Define the path to the VTK file
    base_directory = pwd()
    # Load in simple cubic lattice 41^3
    file = joinpath(base_directory, "post", "centered_simple_cubic_spheres_68921.vtk")
    
    # Check if the file exists
    if !isfile(file)
        error("File not found: $file")
    end

    # Load and process the dataset
    println("Loading file...\n$file")
    data = read_vtk_file(file)
    # execution_time = @elapsed read_vtk_file(file)
    println("Successfully loaded file")
    x_data, y_data, z_data, radii = retrieve_coordinates(data)

    # Convert to cylindrical coordinates
    r_data, theta_data = convert_to_cylindrical(x_data, y_data)

    # Define the lattice size
    a = 40 * 0.004  # Length of one side of the cubic lattice

    # Define Cartesian boundaries
    boundaries_cartesian = Dict(
        :x_min => -a, :x_max => a,
        :y_min => -a, :y_max => a,
        :z_min => -a * 0.2, :z_max => a
    )

    # Define cylindrical boundaries with r_min set to -1 as per your specification
    boundaries_cylindrical = Dict(
        :r_min     => -1,          # Specific meaning in your context
        :r_max     => a,
        :theta_min => π / 2,
        :theta_max => 2π,
        :z_min     => -a * 0.2,
        :z_max     => a
    )

    # Compute packing density in Cartesian coordinates
    packing_density_cartesian = compute_packing_cartesian(
        boundaries=boundaries_cartesian,
        x_data=x_data,
        y_data=y_data,
        z_data=z_data,
        radii=radii
    )

    # Compute packing density in Cylindrical coordinates
    # Measure the execution time using @elapsed
    packing_density_cylindrical, benchmark = benchmarked_call(
        compute_packing_cylindrical
        ;
        boundaries=boundaries_cylindrical,
        r_data=r_data,
        theta_data=theta_data,
        z_data=z_data,
        radii=radii,
        accurate_cylindrical=false
    )


    # Display results
    println("\n=== Packing Density Results ===")
    println("Cartesian Packing Density:    ", packing_density_cartesian)
    println("Cylindrical Packing Density:  ", packing_density_cylindrical)
    println("Simple Cubic Packing Density: ", π / 6)

    println("\n=== Benchmarking Information ===")
    println("Execution Time for compute_packing_cylindrical: $benchmark")
end

function example3()

    # Define the path to the VTK file
    base_directory = pwd()
    file = joinpath(base_directory, "post", "particles_initial.vtk")
    
    # Check if the file exists
    if !isfile(file)
        error("File not found: $file")
    end

    # Load and process the dataset
    println("Loading file...\n$file")
    data = read_vtk_file(file)
    x_data, y_data, z_data, radii = retrieve_coordinates(data)
    r_data, theta_data = convert_to_cylindrical(x_data, y_data)


    # println(size(radii))

    cylinder_radius = 0.03
    cylinder_base_level = 0.0
    cylinder_height = 0.08

    
    ## PACKING DENSITY ##

    boundaries = Dict(
        :r_min => -1,
        :r_max => cylinder_radius*0.9,
        :theta_min => 0,
        :theta_max => pi,
        :z_min => 0.0025,
        :z_max => 0.0325
    )

    packing_density = compute_packing_cylindrical(; boundaries=boundaries,
                                                    r_data=r_data,
                                                    theta_data=theta_data,
                                                    z_data=z_data,
                                                    radii=radii,
                                                    accurate_cylindrical=true)

    # Split data into two types
    types = data[:point_data][:type]

    type_1 = 1.0
    type_2 = 3.0

    radius_inner = nothing

    small_radius = 0.0005
    large_radius = 0.001

    tolerance = 1e-6
    data_1_mask = abs.(radii .- small_radius) .< tolerance
    data_2_mask = abs.(radii .- large_radius) .< tolerance

    # data_1_mask = abs.(types .- type_1) .< tolerance
    # data_2_mask = abs.(types .- type_2) .< tolerance

    data_1 = extract_points(data, data_1_mask)
    data_2 = extract_points(data, data_2_mask)

    num_evals = 100
    target_num_cells_array = LinRange(1, 100000, num_evals)
    SI_vals = zeros(num_evals)

    for (i, target_num_cells) in enumerate(target_num_cells_array)

        SI_vals[i] = calculate_segregation_intensity(data_1,
                                           data_2
                                           ;
                                           cylinder_radius=cylinder_radius,
                                           cylinder_base_level=cylinder_base_level,
                                           cylinder_height=cylinder_height,
                                           target_num_cells=target_num_cells,
                                           output_num_cells=false,
                                           calculate_partial_volumes=true,
                                           clamp_0_to_1=false,
                                           verbose=false)

        print("\rProgress: $(i)/$(num_evals)")
        flush(stdout)

    end

    # Optimise target number of cells by finding the pair of SI values that are closest,
    # and choose the higher SI value along with its corresponding target number.
    if length(SI_vals) < 2
        optimal_target = target_num_cells_array[1]
        optimal_SI = SI_vals[1]
    else
        # Sort the SI values and obtain their original indices
        sorted_idx = sortperm(SI_vals)
        sorted_SI = SI_vals[sorted_idx]
        differences = diff(sorted_SI)
        min_idx = argmin(differences)
        # Choose the higher value of the two adjacent SI values
        chosen_idx = sorted_SI[min_idx] >= sorted_SI[min_idx+1] ? sorted_idx[min_idx] : sorted_idx[min_idx+1]
        optimal_target = target_num_cells_array[chosen_idx]
        optimal_SI = SI_vals[chosen_idx]
    end

    # println("\nOptimal target number of cells: $optimal_target with segregation intensity: $optimal_SI")
    
    # Create the main plot with customised aesthetics.
    p = plot(
        target_num_cells_array, SI_vals,
        xlabel = "Target number of cells",
        ylabel = "Segregation Intensity",
        lw = 2,
        markersize = 3,
        color = :red,
        legend = false,
        grid = true,
        background_color = :white,
        xticks = 0:20000:maximum(target_num_cells_array)
    )

    # # Mark the optimal point on the plot.
    # scatter!(p, [optimal_target], [optimal_SI],
    #         markercolor = :red,
    #         markersize = 8)

    # # Annotate the optimal point.
    # annotate!(p, optimal_target, optimal_SI, text("Optimal", :red, 10))

    # Optionally, adjust the x and y limits for better presentation.
    xlims!(p, 0, 1.1 * maximum(target_num_cells_array))
    ylims!(p, 0, 1)

    # Example values for the ellipse centre and its radii
    centre_x = 53000    # x-coordinate of the ellipse centre
    centre_y = 0.76     # y-coordinate of the ellipse centre
    ellipse_radius_x = 4000   # Horizontal (x-direction) radius of the ellipse
    ellipse_radius_y = 0.1   # Vertical (y-direction) radius of the ellipse

    # Generate points for the ellipse using parametric equations
    t = range(0, stop=2pi, length=200)
    ellipse_x = centre_x .+ ellipse_radius_x .* cos.(t)
    ellipse_y = centre_y .+ ellipse_radius_y .* sin.(t)

    # Add the ellipse to the existing plot with a dashed red line
    plot!(ellipse_x, ellipse_y,
        linestyle = :dash,
        colour = :black,
        lw = 2,
        label = "Highlighted Region")

    # Define the arrow annotation: choose a starting point (for the label) and the tip at the ellipse centre.
    arrow_start = (centre_x + 10000, centre_y - 0.2)  # adjust these offsets as needed
    arrow_tip = (centre_x + 3000, centre_y - 0.06)

    # Annotate the text at the arrow starting location.
    annotate!(arrow_start[1], arrow_start[2], text("Partial volume\ncalculation ceases", 10, :black, :top))

    # Draw a line with an arrow from the annotation point to the ellipse centre.
    plot!([arrow_start[1], arrow_tip[1]], [arrow_start[2], arrow_tip[2]],
        arrow = (:closed, 6),   # the tuple sets the arrowhead style and size
        lw = 2,
        colour = :black)

    # # Annotate the ellipse centre with its coordinates near the centre of the ellipse
    # annotate!(centre_x, centre_y, text("Centre: ($centre_x, $centre_y)", :black, 10))

    # Display the final plot.
    display(p)

    

end

# Call example to test
example1()
