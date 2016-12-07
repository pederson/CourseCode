function [x_hat, Pxx_z] = GaussNewton(dynamics, measurement, t, z, R, xbar, Pxx)
% general Gauss-Newton algorithm
% given a model z = h(x) + epsilon
% where
%       z  -  measurement
%       h(x)  -  model
%       x  -  state vector
%       epsilon - N(0, R) measurement error
%       jacobian(x) - jacobian matrix of h(x) = dh(x)/dx
%
% with some prior - N(xbar, Pxx)



% initialize stuff
xstar = xbar;
dxbar = zeros(size(xbar));
% H = zeros(length(z), length(xbar));
Rinv = inv(R);
Pxxinv = inv(Pxx);

% start looping
for i=1:20
    L = Pxxinv;
    N = L*dxbar;
    
    % add in each measurement
    for k = 1:length(z)
        
        % get the state and STM for this measurement time
        [xst_k, stmk] = dynamics(t(k), xstar);

        % get the measurement info for this state vector
        [hk, Hk] = measurement(xst_k);
        
        % calculate dz
        dz = z(k,:)' - hk;
        
        % update L and N
        L = L + (Hk*stmk)'*Rinv*(Hk*stmk);
        N = N + (Hk*stmk)'*Rinv*dz;

    end
    
    % calculate dx_hat
    dx_hat = L\N;
    
    % calculate new state estimate
    xstar = xstar + dx_hat;
    
    % account for shifting prior
    dxbar = dxbar - dx_hat;
    
end

% calculate Pxx_z
Pxx_z = inv(L);

x_hat = xstar;