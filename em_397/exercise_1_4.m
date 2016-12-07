% Exercise 1.4
% CSE 397 - Bayesian Inverse Problems
% D. Pederson
%

% input parameters
f = @(x) exp(-x);   % pick an analytic function
n = 100;            % number of discrete points
beta = 1;           % std dev of the kernel
sigma = 0.05;       % std dev of noise
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

% Add noise for second part
maxg = max(Yobs);
Yobs = Yobs + sigma*randn(n+1,1);
    

% solve for M
M = A\Yobs;

