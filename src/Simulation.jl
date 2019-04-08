using Compat, Random, Distributions, StatsBase, Gadfly, JuMP, Ipopt, SpecialFunctions, Cairo, Fontconfig

"""
	mom(V::Vector)

Computes the MoM estimator using a vector of gamma-distributed r.v.
"""
function mom(V::Vector)
	μ̂₁ = mean(V)
	μ̂₂ = moment(V, 2, 0)
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
		@variable(model, k, start = (3 - rhs + sqrt((rhs - 3)^2) + 24rhs)/(12rhs))
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

function main()
	Γ =  Gamma(1, 2)
	N = 1000
	M = 5
	x = collect(20:20:N)

	k_mom = zeros(length(x), M)
	s_mom = zeros(length(x), M)
	k_ml = zeros(length(x), M)
	s_ml = zeros(length(x), M)
	index = 1
	for n in x
		for m in 1:M
			V = rand(Γ, n)
			k_mom[index, m], s_mom[index, m] = mom(V)
			k_ml[index, m], s_ml[index, m] = ml(V; mom_init = false)
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

	p_k_mom = plot(x=x, y=k_mom_mean, ymin=k_mom_mean - k_mom_std, ymax=k_mom_mean + k_mom_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt), Guide.xlabel("Size of sample vector"), Guide.ylabel("MoM estimator for k"), Guide.title("Method of moments estimation of the first parameter (k)"));

	p_s_mom = plot(x=x, y=s_mom_mean, ymin=s_mom_mean - s_mom_std, ymax=s_mom_mean + s_mom_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt), Guide.xlabel("Size of sample vector"), Guide.ylabel("MoM estimator for s"), Guide.title("Method of moments estimation of the second parameter (s)"));

	p_k_ml = plot(x=x, y=k_ml_mean, ymin=k_ml_mean - k_ml_std, ymax=k_ml_mean + k_ml_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt), Guide.xlabel("Size of sample vector"), Guide.ylabel("ML estimator for k"), Guide.title("Maximum likelihood estimation of the first parameter (k)"));

	p_s_ml = plot(x=x, y=s_ml_mean, ymin=s_ml_mean - s_ml_std, ymax=s_ml_mean + s_ml_std, Geom.line, Geom.errorbar, Theme(major_label_font="CMU Serif",minor_label_font="CMU Serif",major_label_font_size=16pt,minor_label_font_size=14pt), Guide.xlabel("Size of sample vector"), Guide.ylabel("ML estimator for s"), Guide.title("Maximum likelihood estimation of the second parameter (s)"));

	set_default_plot_size(32cm, 18cm)
	p_k_mom |> PNG("k_mom.png")
	p_s_mom |> PNG("s_mom.png")
	p_k_ml |> PNG("k_ml.png")
	p_s_ml |> PNG("s_ml.png")
end

main()
