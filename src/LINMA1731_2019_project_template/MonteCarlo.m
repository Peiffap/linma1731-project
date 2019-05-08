% sequential_monte_carlo_simple_02PA.m - Started Thu 30 Apr 2015

%   This Matlab script implements a Sequential Monte Carlo (SMC) method for
%   a simple nonlinear dynamical system.
%
%   In the *_incomplete.m version, the goal is to fill out the parts marked by "HIDDEN".
%
%   Reference: Computational Methods in Statistics, Anuj Srivastava, August
%   24, 2009, http://stat.fsu.edu/~anuj/pdf/classes/CompStatII10/BOOK.pdf
%   See in particular problem 3 of section 10.5.


% ** Let us define the dynamical system using the following model:
% x_{t+1} = F(x_t) + Gamma u_t
% y_t = G(x_t) + w_t.
% where x_t in R^m, y_t in R^p,
% x_0 ~ N(mu_x,Sigma_x), u_t ~ N(mu_u,Sigma_u), w_t ~ N(mu_w,Sigma_w)

t_f = 1e2;   % final time. Sugg: 1e2
d_x = 1;  % dimension of state space; must be 1 in this script
d_y = 1;  % dimension of output space; must be 1 in this script
d_u = 1;  % dimension of u; must be 1 in this script
a = .5;  % a (used in F below) should be close to zero for stability
b = 0;  % b (used in F below) should be close to zero, or zero to get a linear system
F = @(x) a*x+b*x^3;  % choice of the function F for the dynamical system
Gamma = 1;
G = @(x) x;  % choice of output function G. Sugg: identity function
mu_x = 0;  % see definition above
Sigma_x = 1;
mu_u = 0;
Sigma_u = 1e-4;
mu_w = 0;
Sigma_w = 1e-0;

sqrt_Sigma_x = sqrt(Sigma_x);
sqrt_Sigma_u = sqrt(Sigma_u);
sqrt_Sigma_w = sqrt(Sigma_w);

out_noise_pdf = @(w) 1/sqrt((2*pi)^d_y*abs(det(Sigma_w))) * exp(-.5*(w-mu_w)'*inv(Sigma_w)*(w-mu_w));  % pdf of the output noise w_t

% ** Simulation: Generate y_t, t=0,..,t_f, that will be used as the
% observations in the SMC algorithm.

x_true = zeros(d_x,t_f +1);  % allocate memory
y_true = zeros(d_y,t_f +1);  % allocate memory

x_true(:,0 +1) = mu_x + sqrt_Sigma_x * randn(d_x,1);  % set true initial state
for t = 0:t_f-1
  u_true = mu_u + sqrt_Sigma_u * randn(d_u,1);  % HIDDEN
  x_true(:,t+1 +1) = F(x_true(:,t +1)) + Gamma * u_true;  % HIDDEN
  w_true = mu_w + sqrt_Sigma_w * randn(d_y,1);  % HIDDEN
  y_true(:,t +1) = G(x_true(:,t +1)) + w_true;  % HIDDEN
end
w_true = mu_w + sqrt_Sigma_w * randn(d_y,1);  % noise on the output at final time t_f
y_true(:,t_f +1) = G(x_true(:,t_f +1)) + w_true;  % output at final time t_f


% *** SEQUENTIAL MONTE CARLO METHOD ***

n = 1e2;   % sample set size. Sugg: 1e2
X = cell(n,t_f +1);   % particles will be stored in X
Xtilde = cell(n,t_f +1);  % to store the predictions

% ** Generate initial sample set {x_0^i,...,x_0^n}:

t = 0;
for i = 1:n
  X{i,t +1} = mu_x + sqrt_Sigma_x * randn(d_x,1);   % we sample from the distribution of x_0  % HIDDEN
end

% ** Start loop on time:

for t = 0:t_f-1
  
  % ** Prediction

  for i = 1:n
    u = mu_u + sqrt_Sigma_u * randn(d_u,1);  % HIDDEN
    Xtilde{i,t+1 +1} = F(X{i,t +1}) + Gamma * u;  % HIDDEN
  end
  
  
  % ** Correction
  
  y = y_true(:,t+1 +1);  % y is the true output at time t+1
  
  weights = zeros(1,n);
  for i=1:n
    weights(i) = out_noise_pdf(y-G(Xtilde{i,t+1 +1}));  % HIDDEN
  end

  % Resample the particles according to the weights:
  if exist('OCTAVE_VERSION') == 0
      % We are using Matlab
      ind_sample = randsample(n,n,true,weights);
  else
      % We are using Octave
      weights_tilde = weights / sum(weights);  % make the weights sum to 1
      [dummy,ind_sample] = histc(rand(n,1),[0 cumsum(weights_tilde)]);
  end

  for i=1:n
    X{i,t+1 +1} = Xtilde{ind_sample(i),t+1 +1};
  end

end  % for t


% ** Visualization
figure(1); clf; hold on
for t = 0:t_f
  % HIDDEN[[  
    % Display particles at each time:  
  for i = 1:n
    plot(t,X{i,t +1},'Marker','o','MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','r')
  end
  % Display true x at each time:
  plot(t,x_true(:,t +1),'kx');
  % Display true y at each time:
  plot(t,y_true(:,t +1),'k>');
  
  % Compute and display sample mean for each time:
  x_mean = zeros(d_x,1);
  for i=1:n
    x_mean = x_mean + X{i,t +1};
  end
  x_mean = x_mean / n;
  plot(t,x_mean,'rx');
  % HIDDEN]]
end
xlabel('t')
ylabel('x_t^i, i=1,...,n')
title('Sequential Monte Carlo experiment')
ax = axis;
axis([-1,t_f,ax(3),ax(4)]);