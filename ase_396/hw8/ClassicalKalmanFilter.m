function [x_hat, Pxx_hat] = ClassicalKalmanFilter(dt, z, R, measurement, ...
                                                  dynamics, control, process_noise, ...
                                                  x, Pxx)
% Classical Kalman Filter update
%
% Given a single measurement z, update the state and covariance

% get the requisite data from the dynamics, control, and noise functions
[~, stmk] = dynamics(dt, x);
[uk, Gk] = control(dt, x);
[vk_bar, Qk, gammak] = process_noise(dt, x);

% time update
xk_bar = stmk*x + Gk*uk + gammak*vk_bar;
Pkbar = stmk*Pxx*stmk' + gammak*Qk*gammak';

% get the requisite data from the measurement function
[~, Hk] = measurement(xk_bar);

% measurement update
Kk = Pkbar*Hk'*inv(Hk*Pkbar*Hk'+R);
x_hat = Kk*(z-Hk*xk_bar) + xk_bar;

Pxx_hat = Pkbar - Kk*Hk*Pkbar;