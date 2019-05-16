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
    d_y = 1;  % dimension of output space; must be 1 in this script
    mu_w = 0;
    t_f = N*ts;
    
    out_noise_pdf = @(w) 1/sqrt((2*pi)^d_y*abs(det(sigma_obs))) * exp(-.5*(w-mu_w)'*inv(sigma_obs)*(w-mu_w));  % pdf of the output noise w_t

    % Initial position guess
    % We compute the sample mean of all the particles we have
    
%     init_pos = zeros(P,2);
%     % y : Px2xN
%     % ye: 1x2xN
%     for i = 1:P
%         init_pos(i,1) = sum(y(i,1,:))/N;
%         init_pos(i,2) = sum(y(i,2,:))/N;
%     end
%     
%     init_posee = zeros(1,2);
%     init_posee(1,1) = sum(ye(1,:,1))/N;
%     init_posee(1,2) = sum(ye(1,:,2))/N;
    
    % Particles will be stores in X
    X  = zeros(P,2,Np,N +1);
    Xe = zeros(1,2,Np,N +1);
    % To store the predictions
    Xtilde  = zeros(P,2,Np,N +1);
    Xtildee = zeros(1,2,Np,N +1);
    
    % How do we get the initial o and oe ?
    % For the moment, let's suppose it is random
    o = y(:, :, 2) - y(:, :, 1) + normrnd(mu_w, sigma_obs, P, 2);  
    norm_o = sqrt(o(:,1).^2 + o(:,2).^2);
    o(:,1) = o(:,1) ./ norm_o;
    o(:,2) = o(:,2) ./ norm_o;
    oe = ye(:, :, 2) - ye(:, :, 1) + normrnd(mu_w, sigma_obs, 1, 2);
    norm_oe = sqrt(oe(:,1).^2 + oe(:,2).^2);
    oe(:,1) = oe(:,1) ./ norm_oe;
    oe(:,2) = oe(:,2) ./ norm_oe;
    
    % How do we now the initial sample set ?
    t = 0;
    for i = 1:Np
        for fish = 1:P
            X(fish,:,i,t +1) = y(fish,:,1) + normrnd(mu_w, sigma_obs, 1, 2);
        end
        Xe(1,:,i,t +1) = ye(1,:,1) + normrnd(mu_w, sigma_obs, 1, 2);
    end
    
    
    % ** Start loop on time:
    for t = 0:N-1
        
        % ** Prediction
        for i = 1:Np
            % Compute the next state for
            [next_x,next_o,next_xe,next_oe] = StateUpdate(X(:,:,i,t +1),o,Xe(1,:,i,t +1),oe,ts,k,s,w);
            o = next_o;
            Xtilde(:,:,i,t+1 +1) = next_x;
            
            oe = next_oe;
            Xtildee(1,:,i, t+1 +1) = next_xe;
        end
        
        % ** Correction

        
        weights  = zeros(P,2,Np);
        weightse = zeros(1,2,Np);
        for i=1:Np
            for fish = 1:P
                weights(fish,1,i) = out_noise_pdf(y(fish,1,t +1)-Xtilde(fish,1,i,t+1 +1));
                weights(fish,2,i) = out_noise_pdf(y(fish,2,t +1)-Xtilde(fish,2,i,t+1 +1));
            end
            weightse(1,1,i) = out_noise_pdf(ye(1,1,t +1)-Xtildee(1,1,i,t+1 +1));
            weightse(1,2,i) = out_noise_pdf(ye(1,2,t +1)-Xtildee(1,2,i,t+1 +1));
        end
        
        ind_sample  = zeros(Np,2,P);
        ind_samplee = zeros(Np,2,1);
        for fish = 1:P
            ind_sample(:,1,fish) = randsample(Np,Np,true,weights(fish,1,:));
            ind_sample(:,2,fish) = randsample(Np,Np,true,weights(fish,2,:));
        end
        ind_samplee(:,1,1) = randsample(Np,Np,true,weightse(1,1,:));
        ind_samplee(:,2,1) = randsample(Np,Np,true,weightse(1,2,:));
        
        for i=1:Np
            for fish = 1:P
                X(fish,1,i,t+1 +1) = Xtilde(fish,1,ind_sample(i),t+1 +1);
                X(fish,2,i,t+1 +1) = Xtilde(fish,2,ind_sample(i),t+1 +1);
            end
            Xe(1,1,i,t+1 +1) = Xtildee(1,1,ind_samplee(i),t+1 +1);
            Xe(1,2,i,t+1 +1) = Xtildee(1,2,ind_samplee(i),t+1 +1);
        end
    end
    
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