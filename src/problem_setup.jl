#=

Run
```julia
include("problem_setup.jl")
```

to generate the same functions

=#
# Circulation
Circulation = OCIM1

# Circulation grid and transport operator
grd, T = Circulation.load()

# Transport operators
T_DIP(p) = T
const DIV = buildDIV(grd)
const Iabove = AIBECS.buildIabove(grd.wet3D, findall(vec(iswet(grd))))
T_POP(p) = transportoperator(grd, w=w(p), DIV=DIV, Iabove=Iabove)
function w(p)
    @unpack w₀, w′ = p
    return @. w₀ + w′ * z
end
iwet = findall(vec(iswet(grd)))
z = ustrip.(grd.depth_3D[iwet])

# Uptake
using WorldOceanAtlasTools, Inpaintings
rawDIPobs3D, nDIPobs3D = WorldOceanAtlasTools.raw_to_grid(grd, 2018, "phosphate", "annual", "1°", "an")
rawDIPobs3D[nDIPobs3D .== 0] .= NaN
DIPobs_3D = fill(NaN, size(grd))
for iz in 1:size(grd)[3]
    DIPobs_3D[:,:,iz] .= inpaint(rawDIPobs3D[:,:,iz])
end
const paintedDIP = DIPobs_3D[iwet]
function U(DIP,p)
    @unpack U₀, α, z₀ = p
    return @. U₀ * DIP^2 / (DIP+α*paintedDIP)^2 * (z≤z₀) * (DIP≥0)
end

# Remin
function R(POP,p)
    @unpack τPOP = p
    return @. POP/τPOP # Warning: Don't allow ∇ₓR = 0 otherwise ∇ₓF is singular at ocean bottom
end

# Net sources and sinks
function G_DIP(DIP, POP, p)
    @unpack xgeo, τgeo = p
    return @. -$U(DIP,p) + $R(POP,p) + (xgeo - DIP) / τgeo
end
function G_POP(DIP, POP, p)
    return @. $U(DIP,p) - $R(POP,p)
end

# Parameters
import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
import AIBECS: @flattenable, flattenable, flatten
@flattenable @units @initial_value struct Params{Tp} <: AbstractParameters{Tp}
    w₀::Tp   | 0.246384 | u"m/d"         | true
    w′::Tp   | 0.134411 | u"m/d/m"       | true
    U₀::Tp   | 1.14631  | u"mmol/m^3/yr" | true
    α::Tp    | 1.0      | u"mol/mol"     | true
    z₀::Tp   | 80.0     | u"m"           | false
    τPOP::Tp | 5.4099   | u"d"           | true
    τgeo::Tp | 1.0      | u"Myr"         | false
    xgeo::Tp | 2.17766  | u"mmol/m^3"    | true
end

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

# Problem setup
nb = sum(iswet(grd))
F, ∇ₓF = state_function_and_Jacobian((T_DIP, T_POP), (G_DIP, G_POP), nb)


