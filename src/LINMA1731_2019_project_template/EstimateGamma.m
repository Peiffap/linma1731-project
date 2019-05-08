function [k,s] = EstimateGamma(noisy_observations)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% The matrix noisy_observaions is 2xn
ts = 0.1; % Time step, should not be written here
[~,n] = size(noisy_observations);

% sqrt( (x2-x1)^2 + (y2-y1)^2 ) / dt
noisy_speed = sqrt((noisy_observations(1,2:n) - noisy_observations(1,1:n-1))^2 + (noisy_observations(2,2:n) - noisy_observations(2,1:n-1))^2)./ts;

% Then here use our function from part 1 to compute k and s

k = 3; s = 2;
end

