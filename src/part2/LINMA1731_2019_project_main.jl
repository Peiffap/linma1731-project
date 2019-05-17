using Compat
using Random, Distributions, StatsBase
using JuMP, Ipopt, SpecialFunctions, LinearAlgebra
using MATLAB
using JLD2

include("EstimateGamma.jl") # Estimation function.
include("ParticleFilter.jl") # Particle filter.

"""
    LINMA1731_2019_project_main

Main script for the LINMA1731 project (2018-2019). The first part of the
script estimates the speed of the fish from noisy measurements and then,
the parameters of the Gamma distribution (using maximum likelihood
estimator) characterising this speed.
Noisy fish trajectories using these parameter estimations are generated.
Finally, a particle filter is applied to track the trajectory of each
fish and of their enemy.

Requirements
------------

  * EstimateGamma: [to implement by yourself using functions from part 1
    of the project!] function that estimates parameters k and s of a Gamma
    distribution from noisy measurements 'noisy_observations'.

  * StateUpdate: function provided by the TA and called by
    GenerateObservations. Useful to update the state vector of the fish and
    the enemy. It implements the state model found in [1].

  * GenerateObservations: function provided by the TA and calling
    recursively StateUpdate and adding Gaussian white noise of variance
    σ_obs^2 to the position and the orientation vectors.

  * ParticleFilter: [to implement by yourself] function tracking the fish
    and the enemy using a particle filter.


Reference
---------

  [1] Huth, A., and Wissel, C. The Simulation of the Movement of Fish
      Schools. Journal of Theoretical Biology 156, 3 (1992), 365�385.

Authors: Charles Wiame and Stephanie Guerit.
         Ported from Matlab to Julia by Louis Navarre and Gilles Peiffer.
Creation: 01-Apr-2019. Last update: 16-May-2019.
Developed: 1.1.0 (2019-01-21)
"""

# Parameters -------------------------------------------------------------------

struct Param
    w::Float64          # Parameter of the size of the FOV.
    P::Int64            # Number of fish.
    N::Int64            # Number of time snapshots.
    Np::Int64           # Number of particles per animal.
    ts::Float64         # Time-step [s].
    sigma_obs::Float64  # Sd of the observation noise
                        # on fish and enemy trajectories.
    k::Float64          # Shape parameter.
    s::Float64          # Scale parameter
end

disp  = true  # Display trajectories.
plot  = false  # Generate arrays needed for plotting (SLOW!)
w     = 20.0
P     = 3
N     = 20
Np    = 20
t_s   = 0.1
σ_obs = 0.2

# Estimation of the parameters of the gamma distribution k and s from noisy
# measurements of the trajectory of one fish -----------------------------------

data = MatFile("noisy_observations.mat") # Open MAT file and return handle.
noisy_observations = get_variable(data, "noisy_observations") # Put in array.
close(data) # Close MAT file.
k̂, ŝ = EstimateGamma(noisy_observations) # TO DEFINE! (keep the same inputs/outputs!)

println(k̂)
println(ŝ)

param = Param(w, P, N, Np, t_s, σ_obs, k̂, ŝ)

# Generate observations --------------------------------------------------------
mat"""
[$x, $xe, $o, $oe, $y, $ye] = GenerateObservations($param);
"""

function test()
    @time mat"""
    ParticleFilter($y, $ye, $param)
    """
end

if disp
    # Particle filtering -------------------------------------------------------
    mat"""
    $x_est, $xe_est = ParticleFilter($y, $ye, $param)
    """

    # Example to display the trajectories. Do not hesitate to adapt it :-)
    mat"""
    for i = 1:$param.N
        cla; hold on
        quiver($x(:,1,i),$x(:,2,i),$o(:,1,i),$o(:,2,i),0,'Marker','o');
        hold on
        quiver($x_est(:,1,i),$x_est(:,2,i),$o(:,1,i),$o(:,2,i),0,'Marker','o', 'color', [1 0 0]);
        hold on
        plot($xe(1,1,i),$xe(1,2,i),'*k');
        hold on
        plot($xe_est(1,1,i),$xe(1,2,i),'*r');
        hold on
        rectangle('Position',[-$param.w -$param.w 2*$param.w 2*$param.w],'EdgeColor','r','LineWidth',3)
        axis([-$param.w-1,$param.w+1,-$param.w-1,$param.w+1]); axis off;
        title(sprintf('time: %3.2f s',i*$param.ts))
        pause(.1);
        hold off;
    end
    """

    MSE_example = (1. /(param.N * param.P) * sum(sum(.√(sum((x - x_est).^2, dims=2)), dims=1), dims=3))[1]
end

# Compute MSE ------------------------------------------------------------------
if plot
    """
        MSE(x, x̂)

    Compute the MSE of a result vector.
    """
    function MSE(x, x̂)
        return (1. /(param.N * param.P) * sum(sum(.√(sum((x - x̂).^2, dims=2)), dims=1), dims=3))[1]
    end

    M = 100  # Number of repetitions.

    Np_vec    = [10 20 50 100]           # Values of Np.
    t_s_vec   = [0.05 0.1 0.2 0.5]       # Values of t_s.
    σ_obs_vec = [0.1 0.2 0.5 1]          # Values of σ_obs.

    MSE_Np    = zeros(length(Np), M)     # MSE for different values of Np.
    MSE_t_s   = zeros(length(t_s), M)    # MSE for different values of t_s.
    MSE_σ_obs = zeros(length(σ_obs), M)  # MSE for different values of σ_obs.

    index = 1

    for Np ∈ Np_vec
        param = Param(w, P, N, Np, t_s, σ_obs, k̂, ŝ)
        for m ∈ 1:M
            # Generate observations ------------------------------------------------
            mat"""
            [$x, $xe, $o, $oe, $y, $ye] = GenerateObservations($param);
            [$x_est, $xe_est] = ParticleFilter($y, $ye, $param);
            """

            MSE_Np[index, m] = MSE(x_est, x̂)
        end
        index += 1
    end

    index = 1

    for t_s ∈ t_s_vec
        param = Param(w, P, N, Np, t_s, σ_obs, k̂, ŝ)
        for m ∈ 1:M
            # Generate observations ----------------------------------------------------
            mat"""
            [$x, $xe, $o, $oe, $y, $ye] = GenerateObservations($param);
            """

            x̂, x̂e = ParticleFilter(y, ye, param)

            MSE_t_s[index, m] = MSE(x, x̂)
        end
        index += 1
    end

    index = 1

    for σ_obs ∈ σ_obs_vec
        param = Param(w, P, N, Np, t_s, σ_obs, k̂, ŝ)
        for m ∈ 1:M
            # Generate observations ----------------------------------------------------
            mat"""
            [$x, $xe, $o, $oe, $y, $ye] = GenerateObservations($param);
            """

            x̂, x̂e = ParticleFilter(y, ye, param)

            MSE_σ_obs[index, m] = MSE(x, x̂)
        end
        index += 1
    end

    @save "MSEArray.jld2" MSE_Np MSE_t_s MSE_σ_obs
end
