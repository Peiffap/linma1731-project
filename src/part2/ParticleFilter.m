function [x_est,xe_est]= ParticleFilter(y,ye,param)

    w = param.w;         % Parameter of the size of the FOV.
    P = param.P;         % Number of fish.
    N = param.N;         % Number of time snapshots.
    Np = param.Np;       % Number of particles per animal.
    ts = param.ts;       % Time-step [s].
    sigma_obs = param.sigma_obs; % Sd of the observation noise
                                 % on fish and enemy trajectories.
    k = param.k;         % Shape parameter.
    s = param.s;         % Scale parameter
    d_y = 1;  % Dimension of output space; must be 1 in this script.
    mu_w = 0;
    t_f = N*ts;

    out_noise_pdf = @(w) 1/sqrt((2*pi)^d_y*abs(det(sigma_obs))) * exp(-.5*(w-mu_w)'*inv(sigma_obs)*(w-mu_w));  % pdf of the output noise w_t

    % Particles will be stored in X(e).
    X  = zeros(P,2,Np,N +1);
    Xe = zeros(1,2,Np,N +1);
    % Predictions are stored in Xtilde(e).
    Xtilde  = zeros(P,2,Np,N +1);
    Xtildee = zeros(1,2,Np,N +1);
    % Orientation for each particle
    o  = zeros(P,2,Np);
    oe = zeros(1,2,Np);

    % The initial orientation is computed as the difference between the first two observations, to which we add some noise.
    % This is a heuristic we found works pretty well.
    for i = 1:Np
        o(:,:,i) = y(:, :, 2) - y(:, :, 1) + normrnd(mu_w, sigma_obs, P, 2);
        norm_o = sqrt(o(:,1).^2 + o(:,2).^2);
        o(:,1,i) = o(:,1) ./ norm_o;
        o(:,2,i) = o(:,2) ./ norm_o;
        oe(1,:,i) = ye(:, :, 2) - ye(:, :, 1) + normrnd(mu_w, sigma_obs, 1, 2);
        norm_oe = sqrt(oe(:,1).^2 + oe(:,2).^2);
        oe(1,1,i) = oe(:,1) ./ norm_oe;
        oe(1,2,i) = oe(:,2) ./ norm_oe;
    end

    % The initial position is guessed by adding a random noise to the observations.
    t = 0;
    for i = 1:Np
        for fish = 1:P
            X(fish,:,i,t +1) = y(fish,:,1) + normrnd(mu_w, sigma_obs, 1, 2);
        end
        Xe(1,:,i,t +1) = ye(1,:,1) + normrnd(mu_w, sigma_obs, 1, 2);
    end


    % Start loop on time; N-2 because of array indexing starting at 1.
    for t = 0:N-2
        % Prediction step.
        % In order to improve performance, we compute the means,
        % and run the simulations with the meaned values.
        for i = 1:Np
        % Compute the next state for
        [next_x,next_o,next_xe,next_oe] = StateUpdate(X(:,:,i,t +1),o(:,:,i),Xe(1,:,i,t +1),oe(:,:,i),ts,k,s,w);
        o(:,:,i) = next_o;
        Xtilde(:,:,i,t+1 +1) = next_x;

        oe(:,:,i) = next_oe;
        Xtildee(1,:,i, t+1 +1) = next_xe;
        end

        % Correction step.
        % Compute the weights for resampling.
        weights  = zeros(P,2,Np);
        weightse = zeros(1,2,Np);
        for i=1:Np
            for fish = 1:P
                weights(fish,1,i) = out_noise_pdf(y(fish,1,t+1 +1)-Xtilde(fish,1,i,t+1 +1));
                weights(fish,2,i) = out_noise_pdf(y(fish,2,t+1 +1)-Xtilde(fish,2,i,t+1 +1));
            end
            weightse(1,1,i) = out_noise_pdf(ye(1,1,t+1 +1)-Xtildee(1,1,i,t+1 +1));
            weightse(1,2,i) = out_noise_pdf(ye(1,2,t+1 +1)-Xtildee(1,2,i,t+1 +1));
        end

        % Sampling.
        % We must take into account the possibility of rounding error on the sample weights leading to a zero vector.
        ind_sample  = zeros(P,Np,2);
        ind_samplee = zeros(1,Np,2);
        for fish = 1:P
            if sum(weights(fish, 1, :)) == 0.0
                ind_sample(fish,:,1) = randsample(Np,Np,true,ones(Np, 1));
            else
                ind_sample(fish,:,1) = randsample(Np,Np,true,weights(fish,1,:));
            end
            if sum(weights(fish, 2, :)) == 0.0
                ind_sample(fish,:,2) = randsample(Np,Np,true,ones(Np, 1));
            else
                ind_sample(fish,:,2) = randsample(Np,Np,true,weights(fish,2,:));
            end
        end
        if sum(weightse(1, 1, :)) == 0.0
            ind_samplee(1,:,1) = randsample(Np,Np,true,ones(Np, 1));
        else
            ind_samplee(1,:,1) = randsample(Np,Np,true,weightse(1,1,:));
        end
        if sum(weightse(1, 2, :)) == 0.0
            ind_samplee(1,:,2) = randsample(Np,Np,true,ones(Np, 1));
        else
            ind_samplee(1,:,2) = randsample(Np,Np,true,weightse(1,2,:));
        end

        % The estimated resampled position is stored in the position vector.
        for i=1:Np
            for fish = 1:P
                X(fish,1,i,t+1 +1) = Xtilde(fish,1,ind_sample(fish,i,1),t+1 +1);
                X(fish,2,i,t+1 +1) = Xtilde(fish,2,ind_sample(fish,i,2),t+1 +1);
            end
            Xe(1,1,i,t+1 +1) = Xtildee(1,1,ind_samplee(1,i,1),t+1 +1);
            Xe(1,2,i,t+1 +1) = Xtildee(1,2,ind_samplee(1,i,2),t+1 +1);
        end
    end

    % Return the mean of the resampled positions as the real estimate.
    x_est = zeros(P,2,N);
    xe_est = zeros(1,2,N);
    for i = 1:N
        for fish = 1:P
            x_est(fish,1,i) = mean(X(fish,1,:,i));
            x_est(fish,2,i) = mean(X(fish,2,:,i));
        end
        xe_est(1,1,i) = mean(Xe(1,1,:,i));
        xe_est(1,2,i) = mean(Xe(1,2,:,i));
    end
end
