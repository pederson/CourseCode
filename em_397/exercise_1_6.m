% Exercise 1.6
% CSE 397 - Bayesian Inverse Problems
% D. Pederson
%

% input parameters
f = @(x) exp(-x);   % pick an analytic function
n = 200;            % number of discrete points
beta = 0.1;           % std dev of the kernel
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
    

% solve using Tikhonov regularization
D = n*(diag(ones(n,1),1) + -1*ones(n+1));
D = D(1:end-1,:);
kap = logspace(-16,2,200);
reg = zeros(size(kap));
for k=1:length(kap)
    Aprime = A'*A + kap(k)*(D')*D;
    bprime = A'*Yobs;

    m = Aprime\bprime;
    reg(k) = norm(Aprime*m-bprime);
end

figure()
loglog(kap, reg, 'ko')
xlabel('\kappa')
ylabel('Misfit |Ax - b|')
title('Minimization of the misfit')

[minreg, imin] = min(reg);
kmin = kap(imin);

hold on
loglog(kmin, minreg, 'ro')


% solution by truncation (Morozov)
xcg = pcg(A, Yobs, sqrt(n+1)*sigma/norm(Yobs), 10000);


% compare to the optimal solution from Tikhonov
Aprime = A'*A+kmin*(D')*D;
bprime = A'*Yobs;
xtik = Aprime\bprime;


figure()
hold on
xvals = (0:n)/n;
plot(xvals, f(xvals), 'k--')
plot(xvals, xtik, 'b*')
plot(xvals, xcg, 'ro')
legend('True Solution','Tikhonov Solution','Morozov solution')


norm(xcg'-f(xvals))
norm(xtik'-f(xvals))