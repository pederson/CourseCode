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

%---------Degenerate prior-------------------

% Construct the L_D matrix
L_D = diag(ones(n+1,1)) - diag(0.5*ones(n,1),1) - diag(0.5*ones(n,1),-1);
L = L_D(2:end-1,:);

% Generate a few standard normal random vectors
nv = 5;
xn = randn(n+1,nv);

% random prior draws (pseudo-inverse)
x = gamma * (L'\xn);

% compute the pointwise variance
Dev = sqrt(gamma^2 * diag(inv(L*L')));

figure
axes('fontsize',12);
plot(s(2:end-1),Dev,'linewidth',2)
legend('Standard deviation','location','best')
hold on
plot(s(2:end-1),x)

%-----------------------------------------------------------------
% Construct the L_D matrix
L_D = diag(ones(n+1,1)) - diag(0.5*ones(n,1),1) - diag(0.5*ones(n,1),-1);

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

%%%-------------Now L_R prior---------------
delta = gamma./ Dev(floor(n/2));
L_R = L_D; 
L_R(1,:) = 0; L_R(1,1) = delta;
L_R(end,:) = 0; L_R(end,end) = delta;


% you should never do this, but we do it anyway for convenience
L_Rinv = inv(L_R);

% random prior draws
x = gamma * (L_Rinv * xn);

% compute the pointwise variance
Dev = sqrt(gamma^2 * diag(L_Rinv * L_Rinv'));

figure
axes('fontsize',12);
plot(s,Dev,'linewidth',2)
legend('Standard deviation','location','best')
hold on
plot(s,x)


%%%-----------Now L_O prior------------------
theta = 0.01;
gamma = 1;
L_N = diag(ones(n+1,1)) - diag(ones(n,1),-1);
m = floor(n/2);
H = ones(n+1,1); H(m) = theta;
L_O = diag(H) * L_N;

% you should never do this, but we do it anyway for convenience
L_Oinv = inv(L_O);

% random prior draws
x = gamma * (L_Oinv * xn);

% compute the pointwise variance
Dev = sqrt(gamma^2 * diag(L_Oinv * L_Oinv'));

figure
axes('fontsize',12);
plot(s,Dev,'linewidth',2)
legend('Standard deviation','location','northwest')
hold on
plot(s,x)
