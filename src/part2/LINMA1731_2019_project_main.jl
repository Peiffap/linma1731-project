using Compat
using Random, Distributions, StatsBase
using JuMP, Ipopt, SpecialFunctions, LinearAlgebra
using PGFPlotsX
using LaTeXStrings
using Distances
using MATLAB

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
         Ported from Matlab to Julia by Louis Navarre and Gilles Peiffer
Creation: 01-Apr-2019. Last update: 09-May-2019.
Developed: 1.1.0 (2019-01-21)
"""

# Parameters --------------------------------------------------------------

struct Param
    w::Float64          # Parameter of the size of the FOV.
    P::Float64          # Number of fish.
    N::Float64          # Number of time snapshots.
    Np::Float64         # Number of particles per animal.
    ts::Float64         # Time-step [s].
    sigma_obs::Float64  # Sd of the observation noise
                        # on fish and enemy trajectories.
    k::Float64          # Shape parameter.
    s::Float64          # Scale parameter
end

disp = true             # Display trajectories.

# Estimation of the parameters of the gamma distribution k and s from noisy
# measurements of the trajectory of one fish ------------------------------

data = MatFile("noisy_observations.mat") # Open MAT file and return handle.
noisy_observations = get_variable(data, "noisy_observations") # Put in array.
close(data) # Close MAT file.

k̂, ŝ = EstimateGamma(noisy_observations) # TO DEFINE! (keep the same inputs/outputs!)

param = Param(20.0, 3.0, 20.0, 20.0, 0.1, 0.2, k̂, ŝ)

# Generate observations ---------------------------------------------------
mat"""
[$x, $xe, $o, $oe, $y, $ye] = GenerateObservations($param);
"""

# Particle filtering ------------------------------------------------------

x_est, xe_est = ParticleFilter(x, xe, param) # TO DEFINE! (keep the same inputs/outputs!)

# Example to display the trajectories. Do not hesitate to adapt it :-)
mat"""
if $disp
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
end
"""

# Compute MSE -------------------------------------------------------------

MSE = (1. /(param.N * param.P) * sum(sum(.√(sum((x - x_est).^2, dims=2)), dims=1), dims=3))[1]
