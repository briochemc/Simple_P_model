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
function U(DIP,p)
    @unpack U₀, k, z₀ = p
    return @. U₀ * DIP/(DIP+k) * (z≤z₀) * (DIP≥0)
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
    w₀::Tp   |  0.64 | u"m/d"         | true
    w′::Tp   |  0.13 | u"m/d/m"       | true
    U₀::Tp   |   1.0 | u"mmol/m^3/yr" | true
    k::Tp    |   0.2 | u"mmol/m^3"    | true
    z₀::Tp   |  80.0 | u"m"           | false
    τPOP::Tp |   5.0 | u"d"           | true
    τgeo::Tp |   1.0 | u"Myr"         | false
    xgeo::Tp |  2.12 | u"mmol/m^3"    | true
end
# Problem setup
nb = sum(iswet(grd))
F, ∇ₓF = state_function_and_Jacobian((T_DIP, T_POP), (G_DIP, G_POP), nb)


