function [xi, w0_m, w0_c, wi] = UnscentedSigmaPoints(xhat, Pxx, alpha, beta, kappa)
%
% Unscented Transform Sigma Points
% D. Pederson
%
% approximate 2n+1 sigma points where n=length(xhat)
%

n = length(xhat);
lam = alpha^2*(n+kappa)-n;
L = chol(Pxx, 'lower');

% form 2n+1 sigma points
xi = zeros(length(xhat), 2*n+1);
xi(:,1) = xhat;
for i=1:n
    xi(:,i+1) = xhat + sqrt(n+lam)*L(:,i);
end
for i=1:n
    xi(:,i+1+n) = xhat - sqrt(n+lam)*L(:,i);
end

% weights
w0_m = lam/(n+lam);
w0_c = (lam/(n+lam) + (1-alpha^2+beta));
wi = 0.5/(n+lam);