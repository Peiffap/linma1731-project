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
    
    out_noise_pdf = @(w) 1/sqrt((2*pi)^d_y*abs(det(sigma_obs))) * exp(-.5*(w-mu_w)'*inv(sigma_obs)*(w-mu_w));  % pdf of the output noise w_t

    % Initial position guess
    % We compute the sample mean of all the particles we have
    
    init_pos = zeros(P,2);
    mm, nn, pp = size(y);
    for i = 1:P
        init_pos(i,1) = sum(y(i,1,:))/pp;
        init_pos(i,2) = sum(y(i,2,:))/pp;
    end
    
    init_posee = zeros(1,2);
    mm, nn = size(y);
    init_posee(1,1) = sum(ye(:,1))/mm; % Ou nn ?
    init_posee(1,2) = sum(ye(:,2))/mm; % Ou nn ?
    
    % Particles will be stores in X
    X = cell(P,2,Np,N +1);
    Xe = cell(2,Np,N +1);
    % To store the predictions
    Xtilde = cell(P,2,Np,N +1);
    Xtildee = cell(2,Np,N +1);
    
    % How do we get the initial o and oe ?
    % For the moment, let's suppose it is random
    o = rand(P,1);
    oe = rand(1,1);
    
    % How do we now the initial sample set ?
    t = 0;
    for i = 1:Np
        for fish = 1:P
            X{fish,:,i,t +1} = x(fish,:,i);
        end
        Xe{:,i,t +1} = xe(:,i);
    end
    
    
    % ** Start loop on time:
    for t = 0:N
        
        % ** Prediction
        for i = 1:Np
            % Compute the next state for i
            [next_x,next_o,next_xe,next_oe] = StateUpdate(X{:,:,i,t +1},o,Xe{:,i,t +1},oe,ts,k,s,w);
            o = next_o;
            Xtilde{:,:,i,t+1 +1} = next_x;
            
            oe = next_oe;
            Xtildee{:,i, t+1 +1} = next_xe;
        end
        
        % ** Correction
        
        yv = y(:,:,t+1 +1); % y is the biased value
        
        weights = zeros(P,2,Np);
        for i=1:Np
            weights(:,2,i) = out_noise_pdf(yv-Xtilde{:,:,i,t+1 +1});
        end
        
        ind_sample = randsample(Np,Np,true, weights);
        
        for i=1:Np
            X{:,:,i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
            Xe{:,i, t+1 +1} = Xtildee{ind_sample(i), t+1 +1};
        end
    end
    
    x_est = zeros(P,2,N);
    xe_est = zeros(2,N);
    for i = 1:N
        for fish = 1:P
            x_est(fish,1,i) = mean(X{fish,1,:,i});
            x_est(fish,2,i) = mean(X{fish,2,:,i});
        end
        xe_est(1,i) = mean(Xe{1,:,i});
        xe_est(2,i) = mean(Xe{2,:,i});
    end

end