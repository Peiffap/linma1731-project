function [next_x,next_o,next_xe,next_oe] = StateUpdate(current_x,current_o,current_xe,current_oe,ts,k_gamma,s_gamma,w)

% StateUpdate is a Matlab function computing the position and the
% orientation of the fish and of the enemy at the next time step.
%
%  StateUpdate(current_x,current_o,current_xe,current_oe,ts,k_gamma,s_gamma,w)
%  computes the position and the orientation of the fish and of the enemy
%  at the next time step based on the current positions and orientations.
%  The state model is based on [1] and adapted to take into account the
%  boundary of the aquarium and the presence of the predator.
%
% Inputs:
%  * current_x:  matrix of size Px2 containing the coordinates of the P
%                fish at current time step.
%  * current_o:  matrix of size Px2 containing the orientations of the P
%                fish at current time step.
%  * current_xe: matrix size 1x2 containing the coordinates of the
%                predator at current time step.
%  * current_oe: matrix size 1x2 containing the orientations of the
%                predator at current time step.
%  * ts:         scalar > 0. Time-step [s].
%  * k_gamma:    scalar > 0. Shape parameter of the Gamma distribution for
%                the fish speed.
%  * s_gamma:    scalar > 0. Rate parameter of the Gamma distribution for
%                the fish speed.
%  * w:          scalar > 0. parameter of the size of the FOV (coordinates
%                of the  corners of the rectangle: (-w,-w),(-w,w),(w,-w),(w,w)).
%
% Outputs:
%  * next_x:   matrix of size Px2 containing the coordinates of the P
%              fish at the next time step.
%  * next_o:   matrix of size Px2 containing the orientations of the P
%              fish at the next time step.
%  * next_xe:  matrix of size 1x2 containing the coordinates of the
%              predator at the time step.
%  * next_oe : matrix of size 1x2 containing the orientation of the
%              predator at the next time step.
%
% References:
%
%  [1] Huth, A., and Wissel, C. The Simulation of the Movement of Fish
%      Schools. Journal of Theoretical Biology 156, 3 (1992), 365–385.
%  [2] Niwa, H.-S. Self-organizing Dynamic Model of Fish Schooling.
%      Journal of Theoretical Biology 171,(1994), 123–136.
%
% Authors: Charles Wiame and Stephanie Guerit.
% Creation: 01-Apr-2019. Last update: 18-Apr-2019.
% Developed: 9.0.0.370719 (R2016a)


% Simulation parameters ---------------------------------------------------
% These parameters values can be used within all your functions.

BL    = 1;            % Body length of the fish
r1    = 0.2*BL;       % Radius of the first area (repulsion area)
r2    = 2.5*BL;       % Radius of the second area (alignement area)
r3    = 5*BL;         % Radius of the third area (attraction area)
r4    = 6*BL;         % Radius of the area in which they see the predator
r5    = 7*BL;         % Radius of the area in which the predator sees the fish
d_m   = 2*BL;         % Security distance to avoid wall collision
ve    = 1*BL;         % Speed of the ennemy
sigma = 7.5;          % Standard deviation of the perturbation on
                      % the turning angles

P = size(current_x,1);             % Number of fish
v = gamrnd(k_gamma,1/s_gamma,P,1); % Fish speed (from a gamma distribution)


% Pre-allocation ----------------------------------------------------------
next_o  = zeros(size(current_o));

% Fish State Update -------------------------------------------------------

% 1. Update the state (position)
next_x  = current_x + current_o.*repmat(v*ts,1,2);
next_x  = min(max(next_x,-w),w);
next_x  = min(max(next_x,-w),w);
next_xe = current_xe + current_oe.*repmat(ve*ts,1,2);

% 2. Update the fish state (orientation)

% Compute the distance matrix D between the fish
D = squareform(pdist(next_x));

for p = 1:P
    
    % Pre-allocation
    beta = zeros(P,1);        % Vector that will contain turning angles
    % between fish j and all other fish
    beta_R = 0;               % Angle with the shark
    danger = 0;               % Boolean flag (true if the fish is close
    % to the predator)
    
    for k = 1:P
        if k ~= p
            
            % Check in which area (repulsion, attraction, etc) is fish
            % k with respect to fish j
            
            if D(p,k) <= r1
                % Repulsion area: fish j goes away from fish k because
                % it is too close.
                a       = computeAngle(current_o(p,:),current_o(k,:));
                tmp     = [a-pi/2 a+pi/2 a-3*pi/2 a+3*pi/2];
                [~,idx] = min(abs(tmp));
                beta(k) = tmp(idx);
                clear a tmp idx
            end
            
            if D(p,k) > r1 && D(p,k) <= r2
                % Alignement area: fish j changes its turning angle
                % to swim parallel to fish k.
                beta(k) = computeAngle(current_o(p,:),current_o(k,:));
            end
            
            if D(p,k) > r2 && D(p,k) <= r3
                % Attraction area: fish j changes its turning angle
                % to swim towards fish k.
                tmp     = (next_x(k,:)-next_x(p,:))...
                    ./norm(next_x(k,:)-next_x(p,:));
                beta(k) = computeAngle(current_o(p,:),tmp);
                clear tmp
            end
            
            if D(p,k) >= r3
                % Searching area: fish j chooses a random angle.
                beta(k) = deg2rad(360*rand(1)-180);
            end
        end
        
    end
    
    if norm(next_x(p,:)-next_xe,2) <= r4
        % If fish j is close to the predactor, it runs away in the
        % opposite direction.
        
        danger = 1;        
        tmp     = -(next_xe-next_x(p,:))./norm(next_xe-next_x(p,:));
        beta_R  = computeAngle(current_o(p,:),tmp);
        clear tmp
    end
    
    % Average decision model: compute the average angle taking into
    % account all the interactions (with the fish and the enemy)
    if danger
        beta_bar = (2*mean(beta)+8*beta_R)/10;
    else
        beta_bar = mean(beta);
    end
    
    % Adding a Gaussian perturbation (due to water flows) to the
    % global turning angle
    A = beta_bar + deg2rad(sigma*randn(1));
    
    
    % Orientation update (using rotation matrix)
    next_o(p,:) = current_o(p,:)*[cos(A) -sin(A) ; sin(A) cos(A)]';
    
    
    % Check if the fish orientation is compatible with the eventual
    % presence of a wall
    
    wall = [false; false; false; false];  % Default value for right, top, left and bottom walls.
    normal_vec = [-1 0; 0 -1; 1 0; 0 1];
    
    if abs(next_x(p,1) - w) <= d_m
        wall(1) = true;
    end
    if abs(next_x(p,2) - w) <= d_m
        wall(2) = true;
    end
    if abs(next_x(p,1) + w) <= d_m
        wall(3) = true;
    end
    if abs(next_x(p,2) + w) <= d_m
        wall(4) = true;
    end
    
    a = pi*ones(length(wall),1);
    for ii = 1:4
        if wall(ii)
            a(ii) = computeAngle(next_o(p,:),-normal_vec(ii,:));
        end
    end
    
    switch length(nonzeros(abs(a) <= round(pi/2,1)))
        case 1
            idxa    = find(abs(a) <= round(pi/2,1));
            tmp     = [a(idxa)-pi/2 a(idxa)+pi/2];
            [~,idx] = min(abs(tmp)); % Choose the angle such that the fish will swim perpendiculary to the wall normal vector.
            A       = tmp(idx);
            
        case 2
            idxa    = find(abs(a) <= round(pi/2,1));
            tmp     = [a(idxa(1))-pi/2 a(idxa(1))+pi/2 a(idxa(2))-pi/2 a(idxa(2))+pi/2];
            [~,idx] = max(abs(tmp)); % Choose the angle such that the fish will swim perpendiculary to the wall normal vector.
            A       = tmp(idx);
            
        otherwise
            A = 0;
    end
    
    % Update the orientation if the fish is too close to a wall
    next_o(p,:) = next_o(p,:)*[cos(A) -sin(A) ; sin(A) cos(A)]';
    
    clear a tmp idx A
    
end

% Predator state update ---------------------------------------------------

Ae_tmp = zeros(P,1);

% Compute the mean angle (determined according to the positions of
% the fish located in the field of view of the enemy).

for p = 1:P
    if norm(next_x(p,:)-next_xe) <= r5 && norm(next_x(p,:)-next_xe) >= 1e-2
        % Enemy goes toward fish
        tmp       = (next_x(p,:)-next_xe)./norm(next_x(p,:)-next_xe);
        Ae_tmp(p) = computeAngle(current_oe,tmp);
        
        clear tmp
    end
end

if isempty(nonzeros(Ae_tmp))
    Ae = deg2rad(360*rand(1)-180);
else
    Ae = mean(nonzeros(Ae_tmp)) + deg2rad(sigma*randn(1));
end

% Orientation update

next_oe = current_oe*[cos(Ae) -sin(Ae) ; sin(Ae) cos(Ae)]';

% Check if the shark orientation is compatible with the eventual
% presence of a wall

wall = [false; false; false; false];  % Default value for right, top, left and bottom walls.
normal_vec = [-1 0; 0 -1; 1 0; 0 1];

if abs(next_xe(1) - w) <= d_m
    wall(1) = true;
end
if abs(next_xe(2) - w) <= d_m
    wall(2) = true;
end
if abs(next_xe(1) + w) <= d_m
    wall(3) = true;
end
if abs(next_xe(2) + w) <= d_m
    wall(4) = true;
end

a = pi*ones(length(wall),1);
for ii = 1:4
    if wall(ii)
        a(ii) = computeAngle(next_oe,-normal_vec(ii,:));
    end
end

switch length(nonzeros(abs(a) <= round(pi/2,1)))
    case 1
        idxa    = find(abs(a) <= round(pi/2,1));
        tmp     = [a(idxa)-pi/2 a(idxa)+pi/2];
        [~,idx] = min(abs(tmp)); % Choose the angle such that the shark will swim perpendiculary to the wall normal vector.
        A       = tmp(idx);
        
    case 2
        idxa = find(abs(a) <= round(pi/2,1));
        tmp  = [a(idxa(1))-pi/2 a(idxa(1))+pi/2 a(idxa(2))-pi/2 a(idxa(2))+pi/2];
        tmp  = sort(tmp,2); % Choose the angle such that the shark will swim perpendiculary to the wall normal vector.
        A    = tmp(end);
        
    otherwise
        A = 0;
end

% Update the orientation if the enemy is too close to a wall
next_oe = next_oe*[cos(A) -sin(A) ; sin(A) cos(A)]';

end

function y = computeAngle(a,b)
% Compute angles in a rigorous way taking the angles into account the four
% quadrants
y = atan2(dot(cross([a 0],[b 0]),[0 0 1]),dot(a,b));
end
