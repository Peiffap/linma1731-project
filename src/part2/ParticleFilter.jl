"""
    ParticleFilter(y, ye, param)

Particle filter for the project.
"""
function ParticleFilter(y, ye, param)
    w     = param.w                   # Parameter of the size of the FOV.
    P     = convert(Int64, param.P)   # Number of fish.
    N     = convert(Int64, param.N)   # Number of time snapshots.
    Np    = convert(Int64, param.Np)  # Number of particles per animal.
    ts    = param.ts                  # Time-step [s].
    σ_obs = param.sigma_obs           # Sd of the observation noise
                                      # on fish and enemy trajectories.
    k     = param.k                   # Shape parameter.
    s     = param.s                   # Scale parameter.

    d_y   = 1  # Dimension of the output space; must be 1 here.
    μ_w   = 0  # Mean of the gaussian error.
    t_f   = N * ts;

    # Probability density function of w_t, evaluated at w.
    α = 1. /√((2π)^d_y * abs(det(σ_obs)))
    out_noise_pdf(w) = α * exp(-0.5(w - μ_w)' * inv(σ_obs) * (w - μ_w))

    # Preallocate storage for the particles.
    X  = zeros(P, 2, Np, N+1)  # Fish particles.
    Xe = zeros(1, 2, Np, N+1)  # Predator particles.

    X̃  = zeros(P, 2, Np, N+1)  # Fish particles -- predictions.
    X̃e = zeros(1, 2, Np, N+1)  # Predator particles -- predictions.

    # Initialize orientations.
    # We assume that the orientation points
    # from the first observation to the second.
    o = y[:, :, 2] - y[:, :, 1]              # Generate random orientations.
    norm_o = .√(o[:, 1].^2 + o[:, 2].^2)     # Compute the norm.
    o[:, 1] = o[:, 1] ./ norm_o              # Normalize.
    o[:, 2] = o[:, 2] ./ norm_o

    oe = ye[:, :, 2] - ye[:, :, 1]           # Generate random orientations.
    norm_oe = .√(oe[:, 1].^2 + oe[:, 2].^2)  # Compute the norm.
    oe[:, 1] = oe[:, 1] ./ norm_oe           # Normalize.
    oe[:, 2] = oe[:, 2] ./ norm_oe

    # Initial sample set is the first observation for each particle.
    t = 0
    for i ∈ 1:Np
        for fish ∈ 1:P
            X[fish, :, i, t+1] = y[fish, :, i]
        end
        Xe[1, :, i, t+1] = ye[1, :, i]
    end

    # Start loop over time.
    println(N-1)
    for t ∈ 0:N-1
    println(t)
        # Prediction step.
        for i ∈ 1:Np
            # Use StateUpdate to get the next position/orientation
            # for the fish and the predator.
            mat"""
            [$next_x, $next_o, $next_xe, $next_oe] = StateUpdate($X(:, :, $i, $t+1), $o, $Xe(1, :, $i, $t+1), $oe, $ts, $k, $s, $w);
            """

            # Next fish orientations and positions.
            o = next_o
            X̃[:, :, i, t+1+1] = next_x

            # Next predator orientations and positions.
            oe = next_oe;
            X̃e[1, :, i, t+1+1] = next_xe
        end

        # Correction step.
        weights  = zeros(P, 2, Np)
        weightse = zeros(P, 2, Np)
        for i ∈ 1:Np
            for fish ∈ 1:P
                # Update weights for each fish.
                weights[fish, 1, i] = out_noise_pdf(y[fish, 1, t+1] - X̃[fish, 1, i, t+1+1])
                weights[fish, 2, i] = out_noise_pdf(y[fish, 2, t+1] - X̃[fish, 2, i, t+1+1])
            end
            # Update weights for the predator.
            weightse[1, 1, i] = out_noise_pdf(ye[1, 1, t+1] - X̃e[1, 1, i, t+1+1])
            weightse[1, 2, i] = out_noise_pdf(ye[1, 2, t+1] - X̃e[1, 2, i, t+1+1])
        end

        ind_sample  = zeros(Np, 2, P)
        ind_samplee = zeros(Np, 2, 1)
        for fish ∈ 1:P
            # Sample for the fish.
            ind_sample[:, 1, fish] = sample(collect(1:Np), pweights(weights[fish, 1, :]), Np)
            ind_sample[:, 2, fish] = sample(collect(1:Np), pweights(weights[fish, 2, :]), Np)
        end
        # Sample for the predator.
        ind_samplee[:, 1, 1] = sample(collect(1:Np), pweights(weightse[1, 1, :]), Np)
        ind_samplee[:, 2, 1] = sample(collect(1:Np), pweights(weightse[1, 2, :]), Np)

        for i ∈ 1:Np
            for fish ∈ 1:P
                # Update real positions with filtered positions for fish.
                X[fish, 1, i, t+1+1] = X̃[fish, 1, convert(Int64, ind_sample[i]), t+1+1]
                X[fish, 2, i, t+1+1] = X̃[fish, 2, convert(Int64, ind_sample[i]), t+1+1]
            end
            # Update real positions with filtered positions for predator.
            Xe[1, 1, i, t+1+1] = X̃e[1, 1, convert(Int64, ind_samplee[i]), t+1+1]
            Xe[1, 2, i, t+1+1] = X̃e[1, 2, convert(Int64, ind_samplee[i]), t+1+1]
        end
    end

    x̂  = zeros(P, 2, N)
    x̂e = zeros(1, 2, N)
    for i ∈ 1:N
        for fish ∈ 1:P
            # Update estimator for the fish.
            x̂[fish, 1, i] = mean(X[fish, 1, :, i])
            x̂[fish, 2, i] = mean(X[fish, 2, :, i])
        end
        # Update estimator for the predator.
        x̂e[1, 1, i] = mean(Xe[1, 1, :, i])
        x̂e[1, 2, i] = mean(Xe[1, 2, :, i])
    end
    return x̂, x̂e
end
