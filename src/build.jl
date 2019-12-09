# This script will generate a simple P model using AIBECS and save
# its steady-state solution in a simple BSON file.

# Run setup script, which should create
# all the parameters and BGC functions,
# as well as the problem
using Pkg
Pkg.activate()
Pkg.instantiate()

using AIBECS, BSON

include("problem_setup.jl")

include("optimization.jl")

# Solve for steady state and save solution
DIP, _ = unpack_tracers(s_optimized, grd)
DIPmismatch = 100sqrt(2AIBECS.mismatch(DIP, μDIPobs, σ²DIPobs, v))

BSON.@save joinpath(@__DIR__, "..", "data", "optimized_solution.bson") s_optimized p_optimized DIPmismatch

#=
To load the optimized model,

using AIBECS, BSON

include("src/problem_setup.jl")

BSON.@load "data/optimized_solution.bson" s_optimized p_optimized

=#
