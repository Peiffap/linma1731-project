function [x_est,xe_est]= ParticleFilter(y,ye,param)

% ParticleFilter is a Matlab function filtering the noisy observations of
% the position of the fish and the predator, using a sequential Monte-Carlo
% filter.
%
%  ParticleFilter(y,ye,param) uses a sequential Monte-Carlo filter to
%  compute an approximation of the true positions of the fish and of the
%  the predator, based on the noisy measurements y (for the fish) and ye
%  (for the predator).
%
% Inputs:
%  * y: matrix of size Px2xN containing the noisy positions of the fish at
%        every time snapshot.
%  * ye: matrix of size 1x2xN containing the noisy positions of the enemy at
%        every time snapshot.
%  * param: Matlab structure containing the following fields:
%    - w:         scalar > 0. parameter of the size of the FOV (coordinates
%                 of the  corners of the rectangle: (-w,-w),(-w,w),(w,-w),(w,w)).
%    - P:         scalar (integer) > 0. Number of fish.
%    - N:         scalar (integer) > 0. Number of time snapshots.
%    - ts:        scalar > 0. Time-step [s].
%    - k:         scalar > 0. Shape parameter of the Gamma distribution for the fish speed.
%    - s:         scalar > 0. Scale parameter of the Gamma distribution for the fish speed.
%    - sigma_obs: scalar > 0. Std of the observation noise on fish and enemy/predator trajectories
%
% Outputs:
%  * x_est:  matrix of size Px2xN containing the approximate positions of the fish at
%        every time snapshot.
%  * xe_est: matrix of size 1x2xN containing the approximate positions of the enemy/predator at
%        every time snapshot.
%
% Requirements:
%  * StateUpdate: function provided by the TA and called by
%    GenerateObservations. Useful to update the state vector of the fish and
%    the enemy. It implements the state model found in [1].
%
% Reference:
%
%  [1] Huth, A., and Wissel, C. The Simulation of the Movement of Fish
%      Schools. Journal of Theoretical Biology 156, 3 (1992), 365--385.
%
% Authors: Louis Navarre and Gilles Peiffer.
% Creation: 07-May-2019. Last update: 17-May-2019.
% Developed: R2017a


% Parameters
w = param.w;         % Parameter of the size of the FOV.
P = param.P;         % Number of fish.
N = param.N;         % Number of time snapshots.
Np = param.Np;       % Number of particles per animal.
ts = param.ts;       % Time-step [s].
sigma_obs = param.sigma_obs; % Sd and mean of the observation noise
                             % on fish and enemy trajectories.
mu_w = 0;                    % Mean of the noise.
k = param.k;         % Shape parameter.
s = param.s;         % Scale parameter
d_y = 1;  % Dimension of output space; must be 1 in this script.

out_noise_pdf = @(w) 1/sqrt((2*pi)^d_y*abs(det(sigma_obs))) * exp(-.5*(w-mu_w)'*inv(sigma_obs)*(w-mu_w));  % Pdf of the output noise w_t.

% Particles will be stored in X(e).
X  = zeros(P,2,Np,N +1);
Xe = zeros(1,2,Np,N +1);
% Predictions are stored in Xtilde(e).
Xtilde  = zeros(P,2,Np,N +1);
Xtildee = zeros(1,2,Np,N +1);
% Orientation for each particle.
o  = zeros(P,2,Np);
oe = zeros(1,2,Np);

% The initial orientation is computed as the difference between the first
% two observations, to which we add some noise to simulate initial
% particles.
% This is a heuristic we found works pretty well.
for i = 1:Np
    % First for the fish.
    % We first compute the difference of the positions between time t=2 and
    % t=1, and we add some noise to it.
    o(:,:,i) = y(:, :, 2) - y(:, :, 1) + normrnd(mu_w, sigma_obs, P, 2);

    % Then, we normalize the orientation vector.
    norm_o = sqrt(o(:,1).^2 + o(:,2).^2);
    o(:,1,i) = o(:,1) ./ norm_o;
    o(:,2,i) = o(:,2) ./ norm_o;

    % Then for the predator.
    % We first compute the difference of the positions between time t=2 and
    % t=1, and we add some noise to it.
    oe(1,:,i) = ye(:, :, 2) - ye(:, :, 1) + normrnd(mu_w, sigma_obs, 1, 2);

    % Then, we normalize the orientation vector.
    norm_oe = sqrt(oe(:,1).^2 + oe(:,2).^2);
    oe(1,1,i) = oe(:,1) ./ norm_oe;
    oe(1,2,i) = oe(:,2) ./ norm_oe;
end % for i = 1:Np

t = 0;

% The only information we have about the initial position of the fish and
% of the predator is the noisy measurement y(e) at time t=1, hence we take
% this knowledge as the initial position.
% To simulate particles from only one measurement, we add some noise to it.
% For each particle.
for i = 1:Np
    % For each fish.
    for fish = 1:P
        X(fish,:,i,t +1) = y(fish,:,1) + normrnd(mu_w, sigma_obs, 1, 2);
    end % for fish = 1:P
    Xe(1,:,i,t +1) = ye(1,:,1) + normrnd(mu_w, sigma_obs, 1, 2);
end % for i = 1:Np

% Iteration over time. Inside this loop, we can see the Monte-Carlo
% sequential method.
% Start loop on time N-2 because of array indexing starting at 1.
for t = 0:N-2

    % ** Prediction step. **

    % For each particle, we call the function StateUpdate.
    % It gives the next position and orientation for each fish and for the predator for the ith particle.
    % We store the results in Xtilde(e), the prediction sets for the fish and the predator, respectively.
    for i = 1:Np
        % Compute the next state.
        [next_x,next_o,next_xe,next_oe] = StateUpdate(X(:,:,i,t +1),o(:,:,i),Xe(1,:,i,t +1),oe(:,:,i),ts,k,s,w);
        % Next orientation for particle i.
        o(:,:,i) = next_o;
        % Prediction sample.
        Xtilde(:,:,i,t+1 +1) = next_x;
        % Next orientation for particle i.
        oe(:,:,i) = next_oe;
        % Prediction sample.
        Xtildee(1,:,i, t+1 +1) = next_xe;
    end % for i = 1:Np

    % ** Correction step. **

    % Compute the weights for resampling.
    % weights: the weights for the fish.
    % weightse: the weights for the predator.
    weights  = zeros(P,2,Np);
    weightse = zeros(1,2,Np);
    for i = 1:Np
        for fish = 1:P
            % For the x-value of the fish.
            weights(fish,1,i) = out_noise_pdf(y(fish,1,t+1 +1)-Xtilde(fish,1,i,t+1 +1));
            % For the y-value of the fish.
            weights(fish,2,i) = out_noise_pdf(y(fish,2,t+1 +1)-Xtilde(fish,2,i,t+1 +1));
        end % for fish = 1:P
        % For the x-value of the predator.
        weightse(1,1,i) = out_noise_pdf(ye(1,1,t+1 +1)-Xtildee(1,1,i,t+1 +1));
        % For the y-value of the predator.
        weightse(1,2,i) = out_noise_pdf(ye(1,2,t+1 +1)-Xtildee(1,2,i,t+1 +1));
    end % for i = 1:Np

    % Sampling.
    % We must take into account the possibility of rounding error on the
    % sample weights leading to a zero vector. If all the weights are
    % rounded to zero, then we call randsample with a vector full of unit
    % values.
    ind_sample  = zeros(P,Np,2);
    ind_samplee = zeros(1,Np,2);
    % For each fish.
    for fish = 1:P
        % For the x-value of the fish.
        if sum(weights(fish, 1, :)) == 0.0
            % Rounding error.
            ind_sample(fish,:,1) = randsample(Np,Np,true,ones(Np, 1));
        else
            ind_sample(fish,:,1) = randsample(Np,Np,true,weights(fish,1,:));
        end % if sum == 0
        % For the y-value of the fish.
        if sum(weights(fish, 2, :)) == 0.0
            % Rounding error.
            ind_sample(fish,:,2) = randsample(Np,Np,true,ones(Np, 1));
        else
            ind_sample(fish,:,2) = randsample(Np,Np,true,weights(fish,2,:));
        end % if sum == 0
    end % for fish = 1:P

    % For the x-value of the predator.
    if sum(weightse(1, 1, :)) == 0.0
        % Rounding error.
        ind_samplee(1,:,1) = randsample(Np,Np,true,ones(Np, 1));
    else
        ind_samplee(1,:,1) = randsample(Np,Np,true,weightse(1,1,:));
    end % if sum == 0
    % For the y-value of the predator.
    if sum(weightse(1, 2, :)) == 0.0
        % Rounding error.
        ind_samplee(1,:,2) = randsample(Np,Np,true,ones(Np, 1));
    else
        ind_samplee(1,:,2) = randsample(Np,Np,true,weightse(1,2,:));
    end % if sum == 0

    % The estimated resampled position is stored in the posterior sets, X(e).
    for i = 1:Np
        for fish = 1:P
            % For the x-value of the fish.
            X(fish,1,i,t+1 +1) = Xtilde(fish,1,ind_sample(fish,i,1),t+1 +1);
            % For the y-value of the fish.
            X(fish,2,i,t+1 +1) = Xtilde(fish,2,ind_sample(fish,i,2),t+1 +1);
        end % for fish = 1:P
        % For the x-value of the predator.
        Xe(1,1,i,t+1 +1) = Xtildee(1,1,ind_samplee(1,i,1),t+1 +1);
        % For the y-value of the predator.
        Xe(1,2,i,t+1 +1) = Xtildee(1,2,ind_samplee(1,i,2),t+1 +1);
    end % for i = 1:Np
end % for t = 0:N-2

% Return the mean of the resampled positions as the real estimate.
x_est = zeros(P,2,N);
xe_est = zeros(1,2,N);
for i = 1:N
    for fish = 1:P
        % For the x-value of the fish
        x_est(fish,1,i) = mean(X(fish,1,:,i));
        % For the y-value of the fish
        x_est(fish,2,i) = mean(X(fish,2,:,i));
    end % for fish = 1:P
    % For the x-value of the predator
    xe_est(1,1,i) = mean(Xe(1,1,:,i));
    % For the y-value of the predator
    xe_est(1,2,i) = mean(Xe(1,2,:,i));
end % for i = 1:N
end % function ParticleFilter
