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

    t_f = N*ts;

    % *** SEQUENTIAL MONTE CARLO FOR THE FISHES ***

    % Currently doing for one fish
    % Should be for all the fishes

    X = cell(N, t_f + 1); % particles wil be stores in X
    Xtiltde = cell(P, N, t_f + 1); % to store the predictions

    % The initial sample set {x_0^i, ..., x_0^n} already exists

    % ** Start loop on time
    
    % !!!!!! We have to compute the initial position of the fishes
    % and of the predator before entering the loop
    % To be able to call the function StateUpdate
    % Guess: initial position is the mean of the observations

    for t = 0:t_f-1
        
        % ** Prediction
        % * Start loop on fishes
        
        for fish = 1:P
            speed = gamrnd(k,s); % The speed is a Gamma distributed RV
            
        end
        
    end

        

end