function [x_est,xe_est]= ParticleFilter(y,ye,param)

w         = param.w;         % Parameter of the size of the FOV
P         = param.P;         % Number of fish
N         = param.N;         % Number of time snapshots
Np        = param.Np;        % Number of particles per animal
ts        = param.ts;        % Time-step [s]
sigma_obs = param.sigma_obs; % Std of the observation noise on fish and 
                             % enemy trajectories

time = ts*N;

% We know the noisy position of the fishes and of the ennemy
% So do we have to compute x_hat ?

end