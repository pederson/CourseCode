function [x_hat, Pxx_hat, dz] = ExtendedKalmanFilter(dt, z, R, measurement, ...
                                                  dynamics, control, process_noise, ...
                                                  x, Pxx)
% Extended Kalman Filter update
%
% Given a single measurement z, update the state and covariance

% get the requisite data from the dynamics, control, and noise functions
[uk, Gk] = control(dt, x);
[vk_bar, Qk, gammak] = process_noise(dt, x);
[xbar_star, stmk] = dynamics(dt, x, uk, vk_bar);


% time update
Pkbar = stmk*Pxx*stmk' + gammak*Qk*gammak';

% get the requisite data from the measurement function
[hk, Hk, dz] = measurement(xbar_star, z);

% measurement update
Kk = Pkbar*Hk'*inv(Hk*Pkbar*Hk'+R);
x_hat = xbar_star + Kk*dz;
Pxx_hat = Pkbar - Kk*Hk*Pkbar;