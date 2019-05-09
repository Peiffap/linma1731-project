


function ParticleFilter(y, ye, param)
    w = param.w         # Parameter of the size of the FOV.
    P = param.P         # Number of fish.
    N = param.N         # Number of time snapshots.
    Np = param.Np       # Number of particles per animal.
    ts = param.ts       # Time-step [s].
    σ = param.sigma_obs # Sd of the observation noise
                        # on fish and enemy trajectories.
    k = param.k         # Shape parameter.
    s = param.s         # Scale parameter
    
    
end