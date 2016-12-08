clear all
% Tan Bui-Thanh, April 2012
% Institute for computational engineering and sciences
% The University of Texas at Austin
% tanbui@ices.utexas.edu

% MCMC algorithms

rand('state',20);
randn('state',18);

% number of samples
N = 1000;

% preallocate the sampling point array
n = 2;
Sample = zeros(n,N);
Samplew = zeros(3,1); % the third one is the weight

% target density
sigma = 0.1; delta = 1;
Pi = @(m) exp(-0.5/sigma^2*(sqrt(sum(m.^2))-1)^2 - 0.5/delta^2* ...
             (m(2)-1)^2);

% proposal variance, and proposal density
gamma = 0.5;
q = @(m,p) 1./sqrt(2*pi*gamma^2)*exp(-0.5/gamma^2*sum((m-p).^2));

% proposal
proposal = @(m) m+gamma*randn(n,1);

disp('starting MCMC...')

% value of pi at the sampling points
Pisampling = zeros(1,N);

% the initial sampling point is at the origin
m = [0,0]';
Samplew(1:n,1) = m; Samplew(n+1,1) = 1;

% target density value of the initial point
PisamplingOld = Pi(m);

sampleCount = 2;

% Metropolis-Hastings loop
while sampleCount <= N,
  % Generate a new point from the current point from the proposal
  % distribution 
  p = proposal(m);
  % Compute the target density
  PisamplingNew = Pi(p);
  % Compute the proposal
  qmp = q(m,p); qpm = q(p,m);
  
  % Compute the acceptance ratio, using log to avoid underflow and
  % overflow
  alpha = log(PisamplingNew*qpm)-log(PisamplingOld* ...
          qmp);
  
  if (alpha > log(rand))
    % accept the new sampling point
    m = p;
    PisamplingOld = PisamplingNew;
    Samplew = [Samplew,ones(n+1,1)];
    Samplew(1:n,end) = m;
  else
    % increase the weight if stay put
    Samplew(n+1,end) = Samplew(n+1,end) + 1;
  end
  
  Sample(:,sampleCount) = m;
  sampleCount = sampleCount+1;
  
end

% average acceptance rate
size(Samplew,2)/N

% estimated mean
mbar = (Samplew(1:n,:) * Samplew(end,:)')/N

figure
axes('fontsize',12);
hold on
plot([-1.5 1.5],[mbar(2), mbar(2)],'k','linewidth',2)
plot([mbar(1) mbar(1)],[-1.5, 1.5],'k','linewidth',2)
%legend('sample mean','location','northeast')

N = 200;
[x,y]=meshgrid(linspace(-1.5,1.5,N));

target = zeros(N,N);
for i = 1:N
  for j = 1:N
    target(i,j) = Pi([x(i,j),y(i,j)]);
  end
end
contour(x,y,target);

N = size(Samplew,2);
t = linspace(0,2*pi,20);
% find the scale
scale = 10*max(Samplew(n+1,:));
% plot a circle at each point with radius proportional to the
% number of staying puts
for i = 1:N
  r = Samplew(n+1,i)/scale;
  x = r*cos(t) + Samplew(1,i);
  y = r*sin(t) + Samplew(2,i);
  fill(x,y,'r')
end
axis equal
axis([-1.5 1.5 -1.5 1.5])

%%%--------- The following two plots only makes sense for large N---

%-------display the trace plot for the first component
figure
axes('fontsize',12);
plot(Sample(1,:))

%--------Compute the correlation length--------------
% subtract out the mean to have zero mean Markov chain
N = size(Sample,2);
Sample0 = Sample - repmat(mean(Sample,2),1,N);
maxLag = 99;
chat = zeros(n,maxLag+1);
for i = 1:n
  c0 = Sample0(i,:)*Sample0(i,:)';
  for k = 0:maxLag
    chat(i,k+1) = sum(Sample0(i,1:N-k) .* Sample0(i,(1:N-k)+k))/c0;
  end
end
figure
axes('fontsize',12);
hold on
odd = 1:2:maxLag+1;
even = 2:2:maxLag+1;
plot(odd,chat(1,odd),'r-o','linewidth',2)
plot(even,chat(2,even),'b-s','linewidth',2)
% $$$ stem(odd,chat(1,odd) , 'filled' , 'r-o');
% $$$ stem(even,chat(2,even) , 'filled' , 'b-s');
legend('1st component','2nd component','location','northeast')
grid  ('on')
xlabel('Lag')
ylabel('Sample Autocorrelation')
title ('Sample Autocorrelation Function (ACF)')
axis([1 100 0 1])
