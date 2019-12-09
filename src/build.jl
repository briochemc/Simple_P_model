# This script will generate a simple P model using AIBECS and save
# its steady-state solution in a simple BSON file.

# Run setup script, which should create
# all the parameters and BGC functions,
# as well as the problem
using AIBECS, BSON

include("problem_setup.jl")

include("optimization.jl")

# Solve for steady state and save solution
BSON.@save joinpath(@__DIR__, "..", "data", "optimized_solution.bson") s_optimized p_optimized

#=
To load the optimized model,

using AIBECS, BSON

include("src/problem_setup.jl")

BSON.@load "data/optimized_solution.bson" s_optimized p_optimized

=#
