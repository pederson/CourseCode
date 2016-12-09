clear all
% Tan Bui-Thanh, April 2012
% Institute for computational engineering and sciences
% The University of Texas at Austin
% tanbui@ices.utexas.edu

% Prior elicitation

rand('state',20);
randn('state',18);

% standard deviation of the innovation
gamma = 0.05;

% Mesh
n = 160;
s = linspace(0,1,n+1)';

%%%%% LAMBDA VALUES
lam = 0.5*ones(n,1); % zero dirichlet boundaries
% lam = 0.5*ones(n,1); lam(1) = 0; lam(end) = 0; % neumann boundaries
%lam = (0:1:n-1)/n;

delta0 = 0.01;
deltan = 0.01;

%-----------------------------------------------------------------
% Construct the L_D matrix
d = ones(n+1, 1); d(1) = delta0;
L_D = gallery('tridiag',-lam, d, lam -1);
L_D(end, end-1) = -1;

% Generate a few standard normal random vectors
nv = 5;
xn = randn(n+1,nv);

% you should never do this, but we do it anyway for convenience
L_Dinv = inv(L_D);

% random prior draws
x = gamma * (L_Dinv * xn);

% compute the pointwise variance
Dev = sqrt(gamma^2 * diag(L_Dinv * L_Dinv'));

figure
axes('fontsize',12);
plot(s,Dev,'linewidth',2)
legend('Standard deviation','location','best')
hold on
plot(s,x)

