function [x_hat, Pxx_hat] = ClassicalKalmanFilter(dynamics, measurement, t, z, R, x, Pxx)
% Classical Kalman Filter update
%
% Given a single measurement z, update the state and covariance

% time update
[xk_bar, stmk] = dynamics(t, x);
Pkbar = stmk*Pxx*stmk';

% measurement update
[hk, Hk] = measurement(x);
Kk = Pkbar*Hk'*inv(Hk*Pkbar*Hk'+R);
x_hat = Kk*(z-Hk*xk_bar) + xk_bar;

Pxx_hat = Pkbar - Kk*Hk*Pkbar;