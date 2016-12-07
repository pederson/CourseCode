function [y_hat, L_hat] = InformationFilter(dt, z, R, measurement, ...
                                                  dynamics, control, process_noise, ...
                                                  y, L)
% Information Filter update
%
% Given a single measurement z, update the state and covariance
% Pxx = inv(L) 
% x = Pxx*y


x = L\y;
% get the requisite data from the dynamics, control, and noise functions
[uk, Gk] = control(dt, x);
[vk_bar, Qk, gammak] = process_noise(dt, x);
[xbar_star, stmk] = dynamics(dt, x, uk, vk_bar);

% time update
Mk = stmk'*L*stmk;
Lk = Mk*gammak*inv(gammak'*Mk*gammak + inv(Qk));
ykbar = (eye(length(y)) - Lk*gammak')*(stmk'*y + Mk*Gk*uk);
Lkbar = Mk - Lk*gammak'*Mk;

% get the requisite data from the measurement function
[hk, Hk] = measurement(xbar_star);

% measurement update
Rinv = inv(R);
L_hat = Hk'*Rinv*Hk + Lkbar;
y_hat = Hk'*Rinv*z + ykbar;
