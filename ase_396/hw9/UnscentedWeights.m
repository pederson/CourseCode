function [w0_m, w0_c, wi] = UnscentedWeights(n, alpha, beta, kappa)
%
% Unscented Transform Weights
% D. Pederson
%

lam = alpha^2*(n+kappa)-n;

% weights
w0_m = lam/(n+lam);
w0_c = (lam/(n+lam) + (1-alpha^2+beta));
wi = 0.5/(n+lam);