% LINMA1731_2019_project_main
%
%  Main script for the LINMA1731 project (2018-2019). The first part of the
%  script estimates the speed of the fish from noisy measurements and then,
%  the parameters of the Gamma distribution (using maximum likelihood
%  estimator) characterising this speed.
%  Noisy fish trajectories using these parameter estimations are generated.
%  Finally, a particle filter is applied to track the trajectory of each
%  fish and of their enemy.
%
% Requirements:
%  * EstimateGamma: [to implement by yourself using functions from part 1 
%    of the project!] function that estimates parameters k and s of a Gamma
%    distribution from noisy measurements 'noisy_observations'. 
%  * StateUpdate: function provided by the TA and called by
%    GenerateObservations. Useful to update the state vector of the fish and
%    the enemy. It implements the state model found in [1].
%  * GenerateObservations: function provided by the TA and calling 
%    recursively StateUpdate and adding Gaussian white noise of variance  
%    sigma_obs^2 to the position and the orientation vectors.
%  * ParticleFilter: [to implement by yourself] function tracking the fish
%    and the enemy using a particle filter.
%
%
% Reference:
%
%  [1] Huth, A., and Wissel, C. The Simulation of the Movement of Fish 
%      Schools. Journal of Theoretical Biology 156, 3 (1992), 365???385.
%
% Authors: Charles Wiame and Stephanie Guerit.
% Creation: 01-Apr-2019. Last update: 18-Apr-2019.
% Developed: 9.0.0.370719 (R2016a)

clear all
clc
close all

% Parameters --------------------------------------------------------------

param.w         = 20;     % Parameter of the size of the FOV
param.P         = 3;     % Number of fish
param.N         = 200;   % Number of time snapshots
param.Np        = 200;    % Number of particles per animal
param.ts        = 0.1;    % Time-step [s]
param.sigma_obs = 0.2;    % Std of the observation noise on fish and 
                          % enemy trajectories
                          
disp = true;              % Display trajectories
                          
% Estimation of the parameters of the gamma distribution k and s from noisy 
% measurements of the trajectory of one fish ------------------------------

load noisy_observations.mat

%[param.k,param.s] = EstimateGamma(noisy_observations); % TO DEFINE! (keep the same inputs/outputs!)
param.k = 3.9540952382110817; param.s = 0.30279168695735004;

% Generate observations ---------------------------------------------------

[x,xe,o,oe,y,ye] = GenerateObservations(param);

tic;
[x_est,xe_est]= ParticleFilter(y,ye,param); % TO DEFINE! (keep the same inputs/outputs!)
toc;

% Exemple to display the trajectories. Do not hesitate to adapt it :-)
if disp
    for i = 1:param.N
        cla; hold on
        quiver(x(:,1,i),x(:,2,i),o(:,1,i),o(:,2,i),0,'Marker','o');
        hold on
        quiver(x_est(:,1,i),x_est(:,2,i),o(:,1,i),o(:,2,i),0,'Marker','o', 'color', [1 0 0]);
        hold on
        plot(xe(1,1,i),xe(1,2,i),'*k');
        hold on
        plot(xe_est(1,1,i),xe(1,2,i),'*r');
        hold on
        rectangle('Position',[-param.w -param.w 2*param.w 2*param.w],'EdgeColor','r','LineWidth',3)
        axis([-param.w-1,param.w+1,-param.w-1,param.w+1]); axis off;
        title(sprintf('time: %3.2f s',i*param.ts))
        pause(.05);
        hold off;
    end
end

% Particle filtering ------------------------------------------------------


% Compute MSE -------------------------------------------------------------

MSE = 1/(param.N*param.P) * sum(sum(sqrt(sum((x-x_est).^2,2)),1),3);