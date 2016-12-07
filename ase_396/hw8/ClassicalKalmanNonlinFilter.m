function [x_hat, Pxx_hat, dx_hat] = ClassicalKalmanNonlinFilter(dt, z, R, measurement, ...
                                                  dynamics, control, process_noise, ...
                                                  x, Pxx, dx)
% Classical Kalman Filter update
%
% Given a single measurement z, update the state and covariance

% get the requisite data from the dynamics, control, and noise functions
[xk_star, stmk] = dynamics(dt, x);
[uk, Gk] = control(dt, x);
[vk_bar, Qk, gammak] = process_noise(dt, x);

% time update
dxk_bar = stmk*dx + Gk*uk + gammak*vk_bar;
Pkbar = stmk*Pxx*stmk' + gammak*Qk*gammak';

% get the requisite data from the measurement function
[hk, Hk] = measurement(xk_star);
dz = z - hk;

% measurement update
Kk = Pkbar*Hk'*inv(Hk*Pkbar*Hk'+R);
dx_hat = Kk*(dz-Hk*dxk_bar) + dxk_bar;

Pxx_hat = Pkbar - Kk*Hk*Pkbar;
x_hat = xk_star + dx_hat;