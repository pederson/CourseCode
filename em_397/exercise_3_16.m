% Exercise 3.16
% CSE 397 - Bayesian Inverse Problems
% D. Pederson
%

% input parameters
f = @(t) 10*(t-0.5).*exp(-0.5*1e2*(t-0.5).^2) -0.8 + 1.6*t;   % pick an analytic function
n = 200;            % number of discrete points
beta = 0.05;           % std dev of the kernel
sigma = 0.01;       % std dev of noise
kappa = 1;          % regularization parameter
%

% the kernel
a = @(s,t,b) 1/sqrt(2*pi*b^2)*exp(-1/(2*b^2)*(t-s)^2);

% construct A matrix
A = zeros(n+1);
for i=0:n
    for j=0:n
        si = i/n;
        sj = j/n;
        
        A(i+1,j+1) = a(si,sj,beta)/n;
    end
end

% construct the noisless observation data
% from the chosen function
Yobs = zeros(n+1,1);
for j=0:n
    for i=0:n
        si = i/n;
        sj = j/n;
        
        Yobs(j+1) = Yobs(j+1) + a(sj,si,beta)*f(si)/n; 
    end
end

% Add noise
maxg = max(Yobs);
Yobs = Yobs + sigma*randn(n+1,1);

% solution by truncation (Morozov)
xcg = pcg(A, Yobs, sqrt(n+1)*sigma/norm(Yobs), 10000);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution by MAP estimation
% Prior flag
PriorFlag = 2; % 1: L_D
               % 2: L_A


% discretize the deblurring kernel
t = (0:1:n)'/(n);
beta = 0.05;
a = 1/sqrt(2*pi*beta^2)*exp(-0.5*(1/beta^2)*t.^2);
A = 1/n*toeplitz(a);

% Truth 
xtrue = 10*(t-0.5).*exp(-0.5*1e2*(t-0.5).^2) -0.8 + 1.6*t;

%%------------additive noise-----------
noise = 100*sigma;       % Noise level in percentages of the max. of noiseless signal
y0 = A*xtrue;    % Noiseless signal
sigma = max(abs(y0))*noise/100;                % STD of the additive noise
y = Yobs;


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
plot(t,xcg,'b')
legend('MAP', 'truth','uncertainty','Morozov estimate','location','best')
plot(t,xmean,'r-','LineWidth',2)
plot(t,xtrue,'k-','LineWidth',1.5)