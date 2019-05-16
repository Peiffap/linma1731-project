using JLD2
using PGFPlotsX
using LaTeXStrings

@load "MSEArray.jld2" MSE_Np MSE_t_s MSE_σ_obs

@pgf p_ratio = Axis(
{
    grid="major",
    title="\\Large Mean square error",
    xlabel=L"\footnotesize Number of particles per fish, \(N_\textnormal{p}\)",
    ylabel=L"\footnotesize Spectral norm of $\cov \hTheta \oslash \mathcal{I}^{-1}(\vartheta) - \left(\begin{smallmatrix}1&1\\1&1\end{smallmatrix}\right)$",
    xmode="log"
},
Plot(Table([:x => N, :y => ratio]))
)
pgfsave("ratio.tex", p_ratio, include_preamble = false)
