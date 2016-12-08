clear all
% Tan Bui-Thanh, April 2012
% Institute for computational engineering and sciences
% The University of Texas at Austin
% tanbui@ices.utexas.edu

% page ranking example

rand('state',20);

% Transition matrix
P = [0 1 0 1/3 0; 
     1/2 0 0 1/3 0; 
     0 0 0 0 1/2; 
     1/2 0 1/2 0 1/2; 
     0 0 1/2 1/3 0];

% number of random variables
n = 5;

I = eye(5);

% number of samples
N = 1500;

x = zeros(n,N);

% initial position
j = 2;
x(:,1) = I(:,j);
pi_j = I(:,j)

for i = 2:N
  % compute the next transtion probability density
  pi_j = P * pi_j;
  
  % compute the CDF
  CDF = cumsum(pi_j);
  
  % draw a random number from the standard uniform distribution
  u = rand(1);
  
  % compute the corresponding random number from PIi by inverse CDF
  k = min(find(CDF>u));
  x(:,i) = I(:,k);
end

% compute the first eigenvector of P corresponding the eigenvalue
% with value 1
[V,D] = eig(P);

% normalize the first eigenvector to become a probability density
a = V(:,1)/sum(V(:,1));

% testing the convergence
norm(a-pi_j)

%testing the convergence of the sample mean.
norm(sum(x,2)/N-a)

%plot the histogram
figure
axes('fontsize',12);
bar([sum(x,2)/N, a])
legend('visiting frequency', 'first eigenvector','location','northeast')