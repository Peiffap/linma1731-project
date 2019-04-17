using Compat
using Random, Distributions, StatsBase
using JuMP, Ipopt, SpecialFunctions, LinearAlgebra
# using Plots; pgfplots()
using PGFPlotsX
using LaTeXStrings

"""
	mom(V::Vector)

Computes the MoM estimator using a vector of gamma-distributed r.v.
"""
function mom(V::Vector)
	μ̂₁ = mean(V) # First moment.
	μ̂₂ = moment(V, 2, 0) # Second moment.

	# Method of moments coefficients.
	k̂ = μ̂₁^2 / (μ̂₂ - μ̂₁^2)
	ŝ = μ̂₂ / μ̂₁ - μ̂₁
	return k̂, ŝ
end

"""
	ml(V::Vector; mom_init = true)

Computes the ML estimator using a vector of gamma-distributed r.v.
"""
function ml(V::Vector; mom_init = false)
	rhs = log(mean(V)) - sum(log(x) for x in V)/length(V) # Right-hand side of the nonlinear equation can be precomputed.

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
	d()

	Plots the output for the answer to question d.
"""
function d()
	Γ =  Gamma(1, 2)
	N = 1000
	M = 500
	x = collect(10:15:N)
	k_mom = zeros(length(x), M)
	s_mom = zeros(length(x), M)
	k_ml = zeros(length(x), M)
	s_ml = zeros(length(x), M)
	index = 1
	for n in x
		for m in 1:M
			V = rand(Γ, n)
			k_mom[index, m], s_mom[index, m] = mom(V) # MoM estimator.
			k_ml[index, m], s_ml[index, m] = ml(V; mom_init = false) # ML estimator.
		end
		index += 1
		println(n)
	end

	k_mom_mean = mean(k_mom[:, i] for i in 1:M)
	s_mom_mean = mean(s_mom[:, i] for i in 1:M)
	k_ml_mean = mean(k_ml[:, i] for i in 1:M)
	s_ml_mean = mean(s_ml[:, i] for i in 1:M)

	k_mom_std = zeros(length(x))
	s_mom_std = zeros(length(x))
	k_ml_std = zeros(length(x))
	s_ml_std = zeros(length(x))
	for i in 1:length(x)
		k_mom_std[i] = std(k_mom[i, :])
		s_mom_std[i] = std(s_mom[i, :])
		k_ml_std[i] = std(k_ml[i, :])
		s_ml_std[i] = std(s_ml[i, :])
	end

	"""
	p_k = plot(
	layer(x=x, y=k_mom_mean, ymin=k_mom_mean - k_mom_std, ymax=k_mom_mean + k_mom_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt, default_color=colorant"blue")),
	layer(x=x, y=k_ml_mean, ymin=k_ml_mean - k_ml_std, ymax=k_ml_mean + k_ml_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt, default_color=colorant"orange")), Guide.xlabel("Size of sample vector"), Guide.ylabel("Mean and standard deviation"), Guide.title("Estimation of the first parameter (k)"), Guide.manual_color_key("Estimators", ["MOM", "ML"], ["blue", "orange"])
	);

	p_s = plot(
	layer(x=x, y=s_mom_mean, ymin=s_mom_mean - s_mom_std, ymax=s_mom_mean + s_mom_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt, default_color=colorant"blue")),
	layer(x=x, y=s_ml_mean, ymin=s_ml_mean - s_ml_std, ymax=s_ml_mean + s_ml_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt, default_color=colorant"orange")), Guide.xlabel("Size of sample vector"), Guide.ylabel("Mean and standard deviation"), Guide.title("Estimation of the second parameter (s)"), Guide.manual_color_key("Estimators", ["MOM", "ML"], ["blue", "orange"])
	);

	p_k |> PNG("k.png")
	p_s |> PNG("s.png")
	"""
end

"""
	f()

	Plots the output for the answer to question f.
"""
function f()
	Γ =  Gamma(1, 2)
	Mf = 10000
	N = [10; 50; 150; 3000]

	k_ml_f = zeros(length(N), Mf)
	s_ml_f = zeros(length(N), Mf)

	index = 1
	for n in N
		for index2 in 1:Mf
			V = rand(Γ, n)
			k_ml_f[index, index2], s_ml_f[index, index2] = ml(V; mom_init=false)
		end
		index += 1
		println(n)
	end

	ratio = zeros(length(N), 1)
	index = 1

	# Values are hardcoded and precomputed for performance reasons.
	ψ = trigamma(1)
	nfisher_inv = [4.0ψ -2.0; -2.0 1.0]/(ψ - 1.0)
	for n in N
		# Compute covariances and store them in a symmetric matrix.
		cov_kk = cov(vec(k_ml_f[index, :]), vec(k_ml_f[index, :]))
		cov_ks = cov(vec(k_ml_f[index, :]), vec(s_ml_f[index, :]))
		cov_ss = cov(vec(s_ml_f[index, :]), vec(s_ml_f[index, :]))
		cov_matrix = [cov_ss cov_ks; cov_ks cov_kk]

		fisher_inverse = nfisher_inv / n

		# Compute the element-wise ratio between
		# the actual covariance and the CRLB
		# and plot its induced 2-norm,
		# and center it around the expected value, that is [1 1; 1 1].
		# This allows us to use the norm to detect convergence
		# (the norm should converge to 0).
		centered_ratio_matrix = (cov_matrix ./ fisher_inverse) .- 1.0
		ratio[index] = opnorm(centered_ratio_matrix, 2) # Spectral norm of the matrix.

		index += 1
	end

	p_ratio = @pgf Axis(
	{
		grid="major",
		title="\\Large Convergence towards the Cramér--Rao bound",
		xlabel="\\footnotesize Size of the sample vector",
		ylabel=L"\footnotesize Spectral norm of $\cov \hTheta \oslash \mathcal{I}^{-1}(\vartheta) - \left(\begin{smallmatrix}1&1\\1&1\end{smallmatrix}\right)$",
		xmode="log"
	},
	Plot(Table([:x => N, :y => ratio]))
	)
	pgfsave("ratio.tex", p_ratio, include_preamble = false)
end

"""
	main()

	Main function to group execution.
"""
function main()
	# d()
	f()
end

main()
