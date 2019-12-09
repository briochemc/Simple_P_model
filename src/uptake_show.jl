using Plots

x = range(0,5,length=100)
xobs = 1
τ = 0.5
k = 0.25

nutrient_restoring(x) = (x - xobs) / τ * (x ≥ xobs)
michaelis_menten_uptake(x) = x / (x + k) * (x ≥ 0)
michaelis_menten_slope(x) = x / τ * x / (x + k) * (x ≥ 0)
mixed_model(x) = x / τ * x / (x + xobs) * (x ≥ 0)
mixed_model_squared(x) = x / τ * x^2 / (x + xobs)^2 * (x ≥ 0)
mixed_model_squared2(x) = x / τ * x^2 / (x + xobs/2)^2 * (x ≥ 0)
mixed_model_exp(x) = x / τ * exp(-xobs / x) * (x ≥ 0)

uptake_functions = [
    nutrient_restoring
    #michaelis_menten_uptake
    michaelis_menten_slope
    mixed_model
    mixed_model_squared
    mixed_model_squared2
    mixed_model_exp
]

plt = plot()

for f in uptake_functions
    plot!(plt, x, f.(x), label=string(f), width=2)
end

plot!(legend=:topleft)

display(plt)

