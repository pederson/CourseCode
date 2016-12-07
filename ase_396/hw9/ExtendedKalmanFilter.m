function [x_hat, Pxx_hat] = ExtendedKalmanFilter(dt, z, R, measurement, ...
                                                  dynamics, control, process_noise, ...
                                                  x, Pxx)
% Extended Kalman Filter update
%
% Given a single measurement z, update the state and covariance

% get the requisite data from the dynamics, control, and noise functions
[xbar_star, stmk] = dynamics(dt, x);
[uk, Gk] = control(dt, x);
[vk_bar, Qk, gammak] = process_noise(dt, x);

% time update
%dxk_bar = stmk*dx + Gk*uk + gammak*vk_bar;
Pkbar = stmk*Pxx*stmk' + gammak*Qk*gammak';

% get the requisite data from the measurement function
[hk, Hk] = measurement(xbar_star);
dz = z - hk;

% measurement update
Kk = Pkbar*Hk'*inv(Hk*Pkbar*Hk'+R);
x_hat = xbar_star + Kk*dz;

Pxx_hat = Pkbar - Kk*Hk*Pkbar;
%x_hat = xk_star + dx_hat;