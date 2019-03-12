using Compat, Random, Distributions, StatsBase, Gadfly

Γ =  Gamma(1, 2)
N = 100000
x = rand(Γ, N)
model = fit(Histogram, x)
plot(model)
