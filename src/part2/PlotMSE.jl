using JLD2
using PGFPlotsX
using LaTeXStrings

@load "MSEArray.jld2" MSE_Np MSE_t_s MSE_σ_obs

@pgf p_ratio = Axis(
{
    grid="major",
    title="\\Large Mean square error",
    xlabel=L"\footnotesize Number of particles per fish, \(N_\textnormal{p}\)",
    ylabel=L"\footnotesize Mean square error $",
    xmode="log"
},
Plot(Table([:x => N, :y => ratio]))
)
pgfsave("ratio.tex", p_ratio, include_preamble = false)
