function [yhat, Pyy] = UnscentedTransform(xhat, Pxx, f, alpha, beta, kappa)
%
% Unscented Transform
% D. Pederson
%
% Nonlinear transformation y = f(x)
% approximated by 2n+1 sigma points where n=length(xhat)
%

n = length(xhat);
lam = alpha^2*(n+kappa)-n;
L = chol(Pxx);

% form 2n+1 sigma points
xi = zeros(length(xhat), 2*n+1);
xi(:,1) = xhat;
for i=1:n
    xi(:,i+1) = xhat + sqrt(n+lam)*L(:,i);
end
for i=1:n
    xi(:,i+1+n) = xhat - sqrt(n+lam)*L(:,i);
end

% propagate the sigma points 
y0 = f(xi(:,1));
yi = zeros(length(y0), 2*n+1);
yi(:,1) = y0;
for i=1:(2*n)
    yi(:,i+1) = f(xi(:,i+1));
end

% get posterior mean
yhat = lam/(n+lam)*yi(:,1);
for i=1:2*n
    yhat = yhat + 0.5/(n+lam)*yi(:,i+1);
end

% get posterior covariance
Pyy = (lam/(n+lam) + (1-alpha^2+beta))*(yi(:,1)-yhat)*(yi(:,1)-yhat)';
for i=1:2*n
    Pyy = Pyy + 0.5/(n+lam)*(yi(:,i+1)-yhat)*(yi(:,i+1)-yhat)';
end