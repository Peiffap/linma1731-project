using JLD2

@load "data/t_s.jld2" MSE_t_s

println(MSE_t_s)
println(maximum(MSE_t_s[1, :]))
