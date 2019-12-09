using Plots

x = range(0,5,length=100)
xobs = 1
τ = 0.5
k = 0.25
α = 1/3

nutrient_restoring(x) = (x - xobs) / τ * (x ≥ xobs)
michaelis_menten_uptake(x) = x / (x + k) * (x ≥ 0)
mixed(x) = x^2 / τ / (x + α * xobs)^2 * (x ≥ 0)

uptake_functions = [
    nutrient_restoring
    michaelis_menten_uptake
    mixed
]

plt = plot()

for f in uptake_functions
    plot!(plt, x, f.(x), label=string(f), width=2)
end

plot!(legend=:topleft)

display(plt)

