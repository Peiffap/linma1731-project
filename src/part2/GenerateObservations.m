function [x,xe,o,oe,y,ye] = GenerateObservations(param)

% GenerateObservations is a Matlab function simulating noisy measurements
% of trajectories.
%
%  GenerateObservations(param) generates true trajectories of fish and
%  their enemy using a slightly modified version of the state model
%  described in [1]. Then, the function adds Gaussian white noise on each
%  coordinate of the position vectors, simulating noisy measurements.
%
% Inputs:
%  * param: Matlab structure containing the following fields:
%    - w:         scalar > 0. parameter of the size of the FOV (coordinates 
%                 of the  corners of the rectangle: (-w,-w),(-w,w),(w,-w),(w,w)).
%    - P:         scalar (integer) > 0. Number of fish.
%    - N:         scalar (integer) > 0. Number of time snapshots.
%    - ts:        scalar > 0. Time-step [s].
%    - k:         scalar > 0. Shape parameter of the Gamma distribution for the fish speed.
%    - s:         scalar > 0. Scale parameter of the Gamma distribution for the fish speed.
%    - sigma_obs: scalar > 0. Std of the observation noise on fish and enemy trajectories
%
% Outputs:
%  * x:  matrix of size Px2xN containing the true positions of the fish at
%        every time snapshot.
%  * xe: matrix of size 1x2xN containing the true positions of the enemy at
%        every time snapshot.
%  * o:  matrix of size Px2xN containing the true orientations of the fish at
%        every time snapshot.
%  * oe: matrix of size 1x2xN containing the true orientations of the enemy at
%        every time snapshot.
%  * y:  matrix of size Px2xN containing the noisy positions of the fish at
%        every time snapshot.
%  * ye: matrix of size 1x2xN containing the noisy positions of the enemy at
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
%      Schools. Journal of Theoretical Biology 156, 3 (1992), 365–385.
%
% Authors: Charles Wiame and Stephanie Guerit.
% Creation: 01-Apr-2019. Last update: 18-Apr-2019.
% Developed: 9.0.0.370719 (R2016a)


% Parameters --------------------------------------------------------------

w         = param.w ;       
P         = param.P;         
N         = param.N;         
ts        = param.ts;        
k         = param.k;         
s         = param.s;         
sigma_obs = param.sigma_obs; 

% Pre-allocation ----------------------------------------------------------

x  = zeros(P,2,N);           % Positions of the fish at each iteration
o  = zeros(P,2,N);           % Orientation (unitary norm) of the fish at each iteration
xe = zeros(1,2,N);           % Positions of the enemy at each iteration
oe = zeros(1,2,N);           % Orientation (unitary norm) of the enemy at each iteration


% Initialisation ----------------------------------------------------------

% Fish
x(:,:,1) = 2*(rand(P,2)-.5)*w/1.5;
o(:,:,1) = rand(P,2) - 0.5;
tmp      = sqrt(sum(o(:,:,1).^2,2));      % Norm of each orientation vector
o(:,1,1) = o(:,1,1)./tmp;
o(:,2,1) = o(:,2,1)./tmp;

% Enemy
xe(:,:,1) = randn(1,2) + 10 ;
oe(:,:,1) = randn(1,2);
oe(:,:,1) = oe(:,:,1)./norm(oe(:,:,1));


% Generate physical trajectories ------------------------------------------

for i = 2:N    
    [x(:,:,i),o(:,:,i),xe(1,:,i),oe(1,:,i)] = ...
        StateUpdate(x(:,:,i-1),o(:,:,i-1),xe(1,:,i-1),oe(1,:,i-1),ts,k,s,w);
end


% Add noise ---------------------------------------------------------------
y  = x  + sigma_obs*randn(P,2,N);
ye = xe + sigma_obs*randn(1,2,N);