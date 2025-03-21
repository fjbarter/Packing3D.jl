using Test
using Packing3D

###########################################################################
# MeshModule Tests
###########################################################################
@testset "MeshModule Tests" begin
    # Test Cartesian Mesh creation with valid parameters
    cart_params = Dict(:x_min => 0.0, :x_max => 1.0, :y_min => 0.0, :y_max => 1.0, :z_min => 0.0, :z_max => 1.0)
    mesh_cart = Mesh(:cartesian; target_num_cells=1000, params=cart_params)
    @test get_total_cells(mesh_cart) > 0
    @test get_cell_volume(mesh_cart) > 0.0
    # Test that get_mesh_bounds returns an array with the expected shape (n_cells x 6)
    bounds = Packing3D.MeshModule.get_mesh_boundaries(mesh_cart)
    @test size(bounds, 2) == 6

    # Test Cylindrical Mesh creation
    cyl_params = Dict(:cylinder_radius => 10.0, :cylinder_base_level => 0.0, :cylinder_height => 20.0)
    cyl_divisions = Dict(:r => 5, :theta => 3, :z => 10)
    mesh_cyl = Mesh(:cylindrical; divisions=cyl_divisions, params=cyl_params)
    @test get_total_cells(mesh_cyl) > 0
    @test get_cell_volume(mesh_cyl) > 0.0

    # Edge Case: Missing parameters should throw an error.
    @test_throws ArgumentError Mesh(:cartesian; target_num_cells=1000, params=Dict(:x_min=>0.0))
end

###########################################################################
# Geometry Tests
###########################################################################
@testset "Geometry Tests" begin
    # Test conversion from Cartesian to Cylindrical
    x_data = [1.0, 0.0, -1.0]
    y_data = [0.0, 1.0, 0.0]
    # Expected: r = [1.0, 1.0, 1.0]; theta values should be approximately 0, π/2, π (modulo 2π)
    r_data, theta_data = convert_to_cylindrical(x_data, y_data)
    @test length(r_data) == length(x_data)
    @test all(isapprox.(r_data, [1.0, 1.0, 1.0], atol=1e-6))
    # For theta, allow modulo differences (e.g. π is equivalent to -π mod 2π)
    function normalize_angle(theta)
      mod(theta, 2*pi)
    end
    @test all(isapprox.(normalize_angle.(theta_data), normalize_angle.([0.0, pi/2, pi]), atol=1e-6))
end

###########################################################################
# IO Tests
###########################################################################
@testset "IO Tests" begin
  # Create a minimal VTK file (ASCII POLYDATA) as a string.
  vtk_string = """
# vtk DataFile Version 5.1
Test VTK file
ASCII
DATASET POLYDATA
POINTS 4 float
0 0 0
1 0 0
1 1 0
0 1 0
"""

  # Create a temporary directory
  temp_dir = mktempdir()
  temp_file = joinpath(temp_dir, "temp.vtk")

  # Write the VTK string to the temporary file.
  open(temp_file, "w") do io
      write(io, vtk_string)
  end

  # Now read from the temporary file.
  data = read_vtk_file(temp_file)
  @test haskey(data, :points)
  @test size(data[:points], 2) == 3
  @test haskey(data, :point_data)

  # Remove the entire temporary directory.
  try
      rm(temp_dir; recursive=true, force=true)
  catch e
      @warn "Failed to remove temporary directory: $temp_dir due to error: $e"
  end
end


###########################################################################
# Data Partitioning Tests (split_data, match_split_data, extract_points)
###########################################################################
@testset "Data Partitioning Tests" begin
    # Create a dummy dataset with 10 points.
    points = [rand(3) for i in 1:10]
    points_matrix = vcat([reshape(p, (1,3)) for p in points]...)
    ids = collect(1:10)
    point_data = Dict(:id => ids, :radius => rand(10))
    data = Dict(:points => points_matrix, :point_data => point_data)

    # Use split_data to partition the data based on x-coordinate.
    data1_ids, data2_ids = split_data(data; split_by=:x)
    @test !isempty(data1_ids)
    @test !isempty(data2_ids)
    # Match the splits to obtain data subsets.
    data1, data2 = match_split_data(data, data1_ids, data2_ids)
    total_points = size(data[:points], 1)
    @test size(data1[:points], 1) + size(data2[:points], 1) == total_points

    # Test extract_points: using a simple mask.
    mask = trues(total_points)
    extracted = extract_points(data, mask)
    @test size(extracted[:points], 1) == total_points
end

###########################################################################
# Public API Functionality Tests
###########################################################################
@testset "Public API Functionality Tests" begin
    # Create dummy particle datasets for two groups.
    N = 50
    x_data = rand(N)
    y_data = rand(N)
    z_data = rand(N)
    radii = 0.01 .+ 0.01 * rand(N)
    # Construct a simple data dictionary.
    data = Dict(:points => hcat(x_data, y_data, z_data),
                 :point_data => Dict(:id => collect(1:N), :radius => radii))
    # Groups 1 and 2: Randomly assigned data to 1 or 2.
    mask_1 = BitVector(rand(Bool, N))
    mask_2 = BitVector(1 .- mask_1)
    data_1 = extract_points(data, mask_1)
    data_2 = extract_points(data, mask_2)
  
    # Define Cartesian boundaries based on data.
    boundaries_cart = [
      minimum(x_data) - 0.01, maximum(x_data) + 0.01,
      minimum(y_data) - 0.01, maximum(y_data) + 0.01,
      minimum(z_data) - 0.01, maximum(z_data) + 0.01
    ]

    params_cart = Dict(
      :x_min => boundaries_cart[1],
      :x_max => boundaries_cart[2],
      :y_max => boundaries_cart[4],
      :y_min => boundaries_cart[3],
      :z_min => boundaries_cart[5],
      :z_max => boundaries_cart[6]
    )

    # Test calculate_packing for Cartesian system.
    packing_cart = calculate_packing(
         ;
         data=data,
         boundaries = boundaries_cart,
         system = :cartesian
    )
    @test isa(packing_cart, Float64)
    @test packing_cart >= 0.0

    # For cylindrical system: convert Cartesian coordinates.
    r_data, theta_data = convert_to_cylindrical(x_data, y_data)
    boundaries_cyl = [
      -1.0, maximum(r_data) + 0.01,
      0.0, 2*pi,
      minimum(z_data) - 0.01, maximum(z_data) + 0.01
    ]

    params_cyl = Dict(
      :cylinder_radius => boundaries_cyl[2],
      :cylinder_base_level => boundaries_cyl[5],
      :cylinder_height => boundaries_cyl[6] - boundaries_cyl[5]
    )

    packing_cyl = calculate_packing(
         ;
         data=data,
         boundaries = boundaries_cyl,
         system = :cylindrical
    )
    @test isa(packing_cyl, Float64)
    @test packing_cyl >= 0.0

    # Test calculate_segregation_intensity for Cartesian system.
    seg_intensity_cart = calculate_segregation_intensity(
         data_1, data_2;
         params = params_cart,
         system = :cartesian,
         target_num_cells = 50
    )
    @test isa(seg_intensity_cart, Float64)
    @test seg_intensity_cart < 1.0

    # Test calculate_segregation_intensity for Cylindrical system.
    seg_intensity_cyl = calculate_segregation_intensity(
         data_1, data_2;
         params = params_cyl,
         system = :cylindrical,
         target_num_cells = 50
    )
    @test isa(seg_intensity_cyl, Float64)
    @test seg_intensity_cyl < 1.0

    # Test calculate_lacey for Cartesian system.
    lacey_cart = calculate_lacey(
         data_1, data_2;
         params = params_cart,
         system = :cartesian,
         target_num_cells = 50
    )
    @test isa(lacey_cart, Float64)
    @test lacey_cart > 0.0

    # Test calculate_lacey for Cylindrical system.
    lacey_cyl = calculate_lacey(
         data_1, data_2;
         params = params_cyl,
         system = :cylindrical,
         target_num_cells = 50
    )
    @test isa(lacey_cyl, Float64)
    @test lacey_cyl > 0.0

    # Test compute_volumes_per_cell.
    # Create a simple Cartesian mesh using the boundaries.
    mesh_cart = Mesh(
         :cartesian;
         target_num_cells = 50,
         params = params_cart
    )

    num_cells = get_total_cells(mesh_cart)

    vols = compute_volumes_per_cell(data_1, data_2; mesh=mesh_cart)
    vols_1, vols_2 = vols
    # println(vols)
    @test length(vols) == 2
    @test length(vols_1) == num_cells
end
