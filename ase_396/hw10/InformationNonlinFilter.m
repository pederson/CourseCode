function [xbar_star, L_hat, dy_hat] = InformationNonlinFilter(dt, z, R, measurement, ...
                                                  dynamics, control, process_noise, ...
                                                  x, L, dy)
% Information Filter update
%
% Given a single measurement z, update the state and covariance
% Pxx = inv(L) 
% x = Pxx*y

% get the requisite data from the dynamics, control, and noise functions
[uk, Gk] = control(dt, x);
[vk_bar, Qk, gammak] = process_noise(dt, x);
[xbar_star, stmk] = dynamics(dt, x, uk, vk_bar);

% time update
Mk = stmk'*L*stmk;
Lk = Mk*gammak*inv(gammak'*Mk*gammak + inv(Qk));
dykbar = (eye(length(x)) - Lk*gammak')*(stmk'*dy + Mk*Gk*uk);
Lkbar = Mk - Lk*gammak'*Mk;

% get the requisite data from the measurement function
[hk, Hk] = measurement(xbar_star);
dz = z - hk;

% measurement update
Rinv = inv(R);
L_hat = Hk'*Rinv*Hk + Lkbar;
dy_hat = Hk'*Rinv*dz + dykbar;

%y_hat = L_hat*xbar_star + dy_hat;
