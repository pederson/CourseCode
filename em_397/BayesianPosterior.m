clear all
% Tan Bui-Thanh, April 2012
% Institute for computational engineering and sciences
% The University of Texas at Austin
% tanbui@ices.utexas.edu

% explore the posterior with smooth priors and without hyper-parameters.

rand('state',20);
randn('state',18);

% Mesh
n = 100; 
s = linspace(0,1,n+1)'; t = s;

% Prior flag
PriorFlag = 2; % 1: L_D
               % 2: L_A


% discretize the deblurring kernel
beta = 0.05;
a = 1/sqrt(2*pi*beta^2)*exp(-0.5*(1/beta^2)*t.^2);
A = 1/n*toeplitz(a);

% Truth 
xtrue = 10*(t-0.5).*exp(-0.5*1e2*(t-0.5).^2) -0.8 + 1.6*t;

%%------------additive noise-----------
noise = 5;       % Noise level in percentages of the max. of noiseless signal
y0 = A*xtrue;    % Noiseless signal
sigma = max(abs(y0))*noise/100;                % STD of the additive noise
y = y0 + sigma*randn(n+1,1);


%%------------Prior construction----------
% standard deviation of the innovation
gamma = 1/n;

% Construct the L_D matrix
if PriorFlag == 1,
  L = diag(ones(n+1,1)) - diag(0.5*ones(n,1),1) - diag(0.5*ones(n,1),-1);
elseif PriorFlag == 2,
  L_D = diag(ones(n+1,1)) - diag(0.5*ones(n,1),1) - diag(0.5*ones(n,1),-1);
  % you should never do this, but we do it anyway for convenience
  L_Dinv = inv(L_D);
  Dev = sqrt(gamma^2 * diag(L_Dinv * L_Dinv'));

  delta = gamma./ Dev(floor(n/2));
  L = L_D; 
  L(1,:) = 0; L(1,1) = delta;
  L(end,:) = 0; L(end,end) = delta;
else
  error('not supported')
end

% Calculating the MAP estimate and posterior variances, by least squares
xmean = [(1/sigma)*A;1/gamma*L]\[(1/sigma)*y;zeros(n+1,1)];
Gamma_post = inv((1/sigma^2)*A'*A + 1/gamma^2*L'*L); 

% Plotting the MAP estimate and the 2*STD envelope

% Defining different shades of blue for plotting
 shades = [176 224 230;
             135 206 235;
             135 206 255;
             126 192 238;
             108 166 205];
   shades = 1/255*shades;

STD = sqrt(diag(Gamma_post));
xhigh = xmean + 2*STD;
xlow = xmean - 2*STD;

figure
axes('fontsize',12);
plot(t,xmean,'r-','LineWidth',2), hold on
plot(t,xtrue,'k-','LineWidth',1.5)
fill([t;t(n+1:-1:1)],[xlow;xhigh(n+1:-1:1)],shades(1,:))
legend('MAP', 'truth','uncertainty','location','best')
plot(t,xmean,'r-','LineWidth',2)
plot(t,xtrue,'k-','LineWidth',1.5)