using Compat
using Random, Distributions, StatsBase
using JuMP, Ipopt, SpecialFunctions, LinearAlgebra

"""
	mom(V::Array{Float64, 1})

Computes the MoM estimator using a vector of gamma-distributed r.v.
"""
function mom(V::Array{Float64, 1})
	μ̂₁ = mean(V) # First moment.
	μ̂₂ = moment(V, 2, 0) # Second moment.

	# Method of moments coefficients.
	k̂ = μ̂₁^2 / (μ̂₂ - μ̂₁^2)
	ŝ = μ̂₂ / μ̂₁ - μ̂₁
	return k̂, ŝ
end

"""
	ml(V::Array{Float64, 1}; mom_init = true)

Computes the ML estimator using a vector of gamma-distributed r.v.
"""
function ml(V::Array{Float64, 1}; mom_init = false)
	rhs = log(mean(V)) - sum(log.(V))/length(V) # Right-hand side of the nonlinear equation can be precomputed.

	# The nonlinear optimisation approach uses the fact that
	# if one optimises a constant function,
	# using the nonlinear equation as a constraint,
	# then the result will be the solution of the equation.
	model = Model(with_optimizer(Ipopt.Optimizer, print_level = 0))
	if mom_init # Initial guess using method of moments estimation.
		@variable(model, k, start = mom(V)[1])
	else
		# Minka, Thomas P. (2002). "Estimating a Gamma distribution",
		# (https://tminka.github.io/papers/minka-gamma.pdf).
		@variable(model, k, start = (3.0 - rhs + sqrt((rhs - 3.0)^2) + 24rhs)/(12rhs))
	end
	@NLobjective(model, Min, 1.0)

	# Adding a "user-defined" function in constraints requires some magic...
	my_digamma(k) = digamma(k)
	my_trigamma(k) = trigamma(k)
	my_tetragamma(k) = polygamma(2, k)
	register(model, :my_digamma, 1, my_digamma, my_trigamma, my_tetragamma)
	@NLconstraint(model, log(k) - my_digamma(k) == rhs)
	JuMP.optimize!(model)

	k̂ = JuMP.value(k)
	ŝ = mean(V)/k̂
	return k̂, ŝ
end

"""
	EstimateGamma(obs::Array{Float64, 2})

Computes the maximum likelihood estimators of the distribution of speeds
while taking as argument the vector of 2D positions.
"""
function EstimateGamma(obs::Array{Float64, 2})
	ts = 0.1 # Time step [s].
	_, n = size(obs)

	v = .√((obs[1, 2:n] - obs[1, 1:n-1]).^2 + (obs[2, 2:n] - obs[2, 1:n-1]).^2) / ts # Compute the speeds as explained in the report.

	return ml(v)
end
