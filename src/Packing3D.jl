# Packing3D.jl # CASE SENSITIVE

# This is the parent module for the packing3d package
module Packing3D

# Include each module
include("geometry.jl")
include("io.jl")
include("mesh.jl")
include("utils.jl")
include("cartesian.jl")
include("cylindrical.jl")

# Bring in submodules (local imports)
using .Geometry
using .IO
using .MeshModule
using .Utils
using .Cartesian
using .Cylindrical

# Exports already dealt with inside modules, repeated here for clarity
export Mesh,
       compute_packing_cartesian,
       compute_packing_cylindrical,
       calculate_segregation_intensity,
       calculate_lacey,
       convert_to_cylindrical,
       read_vtk_file,
       retrieve_coordinates,
       extract_points


end # module Packing3D