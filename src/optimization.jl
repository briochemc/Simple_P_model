# 1st solve (needed?)
nt = 2
p = Params()
@unpack xgeo = p
x = xgeo * ones(nt * nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)
s = solve(prob, CTKAlg(), preprint="1st solve ").u

# optimize model

# Assign a lognormal prior to each based on initial value
using Distributions
import AIBECS: @prior, prior
function prior(::Type{T}, s::Symbol) where {T<:Params}
    if flattenable(T, s)
        μ = log(ustrip(upreferred(initial_value(T, s) * units(T, s))))
        return LogNormal(μ ,1.0)
    else
        return nothing
    end
end
prior(::T, s::Symbol) where {T<:AbstractParameters} = prior(T,s)
prior(::Type{T}) where {T<:AbstractParameters} = Tuple(prior(T,s) for s in AIBECS.symbols(T))
prior(::T) where {T<:AbstractParameters} = prior(T)

# WOA data
using WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, 2018, "phosphate", "annual", "1°", "an")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
μx = (μDIPobs, missing)
σ²x = (σ²DIPobs, missing)

# param weights
ωs = [1.0, 0.0] # the weight for the mismatch (weight of POP = 0)
ωp = 1e-4

# objective functions
v = ustrip.(vector_of_volumes(grd))
f, ∇ₓf, ∇ₚf = generate_objective_and_derivatives(ωs, μx, σ²x, v, ωp)

# Use F1 for gradient and Hessian
using F1Method
mem = F1Method.initialize_mem(s, p)
objective(p) = F1Method.objective(f, F, ∇ₓF, mem, p, CTKAlg(), preprint="obj ")
gradient(p) = F1Method.gradient(f, F, ∇ₓf, ∇ₓF, mem, p, CTKAlg(), preprint="grad ")
hessian(p) = F1Method.hessian(f, F, ∇ₓf, ∇ₓF, mem, p, CTKAlg(), preprint="hess ")

# change of variables
λ2p = subfun(typeof(p))
∇λ2p = ∇subfun(typeof(p))
∇²λ2p = ∇²subfun(typeof(p))
p2λ = invsubfun(typeof(p))
λ = p2λ(p)

#===========================================================================
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
# EDIT HERE
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
===========================================================================#

# variable-changed objectiv function
function obj(λ)
    show(λ2p(λ))
    return objective(λ2p(λ))
end
using LinearAlgebra
function grad(λ)
    show(λ2p(λ))
    return gradient(λ2p(λ)) * Diagonal(∇λ2p(λ))
end
function hess(λ)
    show(λ2p(λ))
    ∇p = Diagonal(∇λ2p(λ)) # for variable change
    ∇²p = Diagonal(∇²λ2p(λ)) # for variable change
    G = vec(gradient(λ2p(λ)))
    H = hessian(λ2p(λ))
    return ∇p * H * ∇p + Diagonal(G) * ∇²p
end

using Optim
m = length(p)
grad(s, λ) = s[1:m] .= vec(grad(λ))
hess(s, λ) = s[1:m,1:m] .= hess(λ)
opt = Optim.Options(store_trace = false, show_trace = true, extended_trace = false)
results = optimize(obj, grad, hess, λ, NewtonTrustRegion(), opt)

p_optimized = λ2p(results.minimizer)
prob_optimized = SteadyStateProblem(F, ∇ₓF, s, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg()).u
