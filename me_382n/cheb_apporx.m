function [values, ak] = cheb_apporx ( function_evaluator_handle, points, M)
% construct L2 approximation of a function using
% chebyshev basis functions
% The function is defined on -1 to 1 for x which means
% that it's 0 to pi on s
%N = length(points);

% define coordinate transformation
s = @(x) acos(-x);
q = @(s) function_evaluator_handle(-cos(s));

% define the basis functions
bks = @(s, k) cos(k.*s);
bkx = @(x, k) cos(k*s(x));

si = [0:M]'*pi/M;
ki = [0:M]';

% evaluate coefficients by the trapezoid rule
ak = arrayfun(@(k) 2/M*(-(q(0) + q(pi)*(-1).^k)/2 + sum(q(si).*bks(si, k))), ki);


% compute the approximation of the values
values = arrayfun(@(p) ak(1)/2 + dot(ak(2:end),bkx(p, ki(2:end))), points);




end
