using JLD2
using PGFPlotsX
using LaTeXStrings

@load "data/np.jld2" MSE_Np Np_vec
@load "data/np_y.jld2" MSE_Np_y
@load "data/t_s.jld2" MSE_t_s t_s_vec
@load "data/t_s_y.jld2" MSE_t_s_y
@load "data/sigma_obs.jld2" MSE_σ_obs σ_obs_vec
@load "data/sigma_obs_y.jld2" MSE_σ_obs_y

M = 15                         # Number of repetitions.
Np_len = length(Np_vec)        # Number of tested values for Np.
t_s_len = length(t_s_vec)      # Number of tested values for t_s.
σ_obs_len = length(σ_obs_vec)  # Number of tested values for σ_obs.

Np_mean      = mean(MSE_Np[:, i] for i ∈ 1:M)
Np_y_mean    = mean(MSE_Np_y[:, i] for i ∈ 1:M)
t_s_mean     = mean(MSE_t_s[:, i] for i ∈ 1:M)
t_s_y_mean   = mean(MSE_t_s_y[:, i] for i ∈ 1:M)
σ_obs_mean   = mean(MSE_σ_obs[:, i] for i ∈ 1:M)
σ_obs_y_mean = mean(MSE_σ_obs_y[:, i] for i ∈ 1:M)

Np_std       = zeros(Np_len)
Np_y_std     = zeros(Np_len)
t_s_std      = zeros(t_s_len)
t_s_y_std    = zeros(t_s_len)
σ_obs_std    = zeros(σ_obs_len)
σ_obs_y_std  = zeros(σ_obs_len)

for i ∈ 1:Np_len
    Np_std[i]      = std(MSE_Np[i, :])
    Np_y_std[i]    = std(MSE_Np_y[i, :])
end

for i ∈ 1:t_s_len
    t_s_std[i]     = std(MSE_t_s[i, :])
    t_s_y_std[i]   = std(MSE_t_s_y[i, :])
end

for i ∈ 1:σ_obs_len
    σ_obs_std[i]   = std(MSE_σ_obs[i, :])
    σ_obs_y_std[i] = std(MSE_σ_obs_y[i, :])
end

@pgf p_mse_np = Axis(
{
    grid="major",
    title=L"\Large Influence of $N_\textnormal{p}$ on the mean squared error",
    xlabel=L"\footnotesize Number of particles per fish, $N_\textnormal{p}$",
    ylabel="\\footnotesize Mean squared error",
    xmode="log",
    "legend pos=north east",
},
PlotInc({"blue, error bars/.cd, y dir=both,y explicit,
error bar style={line width=1pt},
error mark options={
    rotate=90,
    blue,
    mark size=4pt,
    line width=1pt
}"}, Table({"y error=error"}, [:x => Np_vec', :y => Np_mean, :error => Np_std])),
PlotInc({"red, error bars/.cd, y dir=both,y explicit,
error bar style={line width=1pt},
error mark options={
    rotate=90,
    red,
    mark size=4pt,
    line width=1pt
}"}, Table({"y error=error"}, [:x => Np_vec', :y => Np_y_mean, :error => Np_y_std])),
Legend([L"$E_\textnormal{MSE}\big(\hat{\bm{x}}\big)$", L"$E_\textnormal{MSE}\big(\bm{y}\big)$"])
)
pgfsave("../../report/img/mse_np.tex", p_mse_np, include_preamble = false)

@pgf p_mse_t_s = Axis(
{
    grid="major",
    title=L"\Large Influence of $t_\textnormal{s}$ on the mean squared error",
    xlabel=L"\footnotesize Sampling period, $t_\textnormal{s}$",
    ylabel="\\footnotesize Mean squared error",
    xmode="log",
    "legend pos=north west",
},
PlotInc({"blue, error bars/.cd, y dir=both,y explicit,
error bar style={line width=1pt},
error mark options={
    rotate=90,
    blue,
    mark size=4pt,
    line width=1pt
}"}, Table({"y error=error"}, [:x => t_s_vec', :y => t_s_mean, :error => t_s_std])),
PlotInc({"red, error bars/.cd, y dir=both,y explicit,
error bar style={line width=1pt},
error mark options={
    rotate=90,
    red,
    mark size=4pt,
    line width=1pt
}"}, Table({"y error=error"}, [:x => t_s_vec', :y => t_s_y_mean, :error => t_s_y_std])),
Legend([L"$E_\textnormal{MSE}\big(\hat{\bm{x}}\big)$", L"$E_\textnormal{MSE}\big(\bm{y}\big)$"])
)
pgfsave("../../report/img/mse_t_s.tex", p_mse_t_s, include_preamble = false)

@pgf p_mse_σ_obs = Axis(
{
    grid="major",
    title=L"\Large Influence of $\sigma_\textnormal{obs}$ on the mean squared error",
    xlabel=L"\footnotesize Standard deviation of the noise on the observations, $\sigma_\textnormal{obs}$",
    ylabel="\\footnotesize Mean squared error",
    xmode="log",
    "legend pos=north west",
},
PlotInc({"blue, error bars/.cd, y dir=both,y explicit,
error bar style={line width=1pt},
error mark options={
    rotate=90,
    blue,
    mark size=4pt,
    line width=1pt
}"}, Table({"y error=error"}, [:x => σ_obs_vec', :y => σ_obs_mean, :error => σ_obs_std])),
PlotInc({"red, error bars/.cd, y dir=both,y explicit,
error bar style={line width=1pt},
error mark options={
    rotate=90,
    red,
    mark size=4pt,
    line width=1pt
}"}, Table({"y error=error"}, [:x => σ_obs_vec', :y => σ_obs_y_mean, :error => σ_obs_y_std])),
Legend([L"$E_\textnormal{MSE}\big(\hat{\bm{x}}\big)$", L"$E_\textnormal{MSE}\big(\bm{y}\big)$"])
)
pgfsave("../../report/img/mse_sigma_obs.tex", p_mse_σ_obs, include_preamble = false)
