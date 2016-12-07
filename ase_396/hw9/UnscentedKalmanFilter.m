function [x_hat, Pxx_hat] = UnscentedKalmanFilter(dt, z, R, measurement, ...
                                                  dynamics, control, process_noise, ...
                                                  x, Pxx)
% Unscented Kalman Filter update
%
% Given a single measurement z, update the state and covariance

[v,Q] = process_noise(dt, x);
u = control(dt, x);

n = length(x);
d = length(v);
m = length(z);

% generate the sigma points
P = zeros(n+d+m, n+d+m);
P(1:n,1:n) = Pxx;
P((n+1):(n+d), (n+1):(n+d)) = Q;
P((n+d+1):end, (n+d+1):end) = R;

%{
P
Pxx
Q
R
%}

[xi, w0_m, w0_c, wi] = UnscentedSigmaPoints([x; v; zeros(length(z),1)],P, 1, 3, 3-(n+d+m));
[nstate, nsig] = size(xi);

%{
n
d
m
nstate 
nsig
%}

                                          
% do time update from dynamics
for i=1:nsig
    xi(1:n,i) = dynamics(dt, xi(1:n,i), u, xi((n+1):(n+d),i));
end

% calculate xbar and Pxxbar
xbar = w0_m*xi(1:n,1);
for i=2:nsig
    xbar = xbar + wi*xi(1:n,i);
end
Pxxbar = w0_c*(xi(1:n,1)-xbar)*(xi(1:n,1)-xbar)';
for i=2:nsig
    Pxxbar = Pxxbar + wi*(xi(1:n,i)-xbar)*(xi(1:n,i)-xbar)';
end

% use measurement model to generate Zi
Zi = zeros(length(z),nsig);
for i=1:nsig
    Zi(:,i) = measurement(xi(1:n,i)) + xi(n+d+1:end,i);
end

% calculate zbar and Pzzbar
zbar = w0_m*Zi(:,1);
for i=2:nsig
    zbar = zbar + wi*Zi(:,i);
end
Pzzbar = w0_c*(Zi(:,1)-zbar)*(Zi(:,1)-zbar)';
for i=2:nsig
    Pzzbar = Pzzbar + wi*(Zi(:,i)-zbar)*(Zi(:,i)-zbar)';
end

% calculate cross covariance Pxzbar
Pxzbar = w0_c*(xi(1:n,1)-xbar)*(Zi(:,1)-zbar)';
for i=2:nsig
    Pxzbar = Pxzbar + wi*(xi(1:n,i)-xbar)*(Zi(:,i)-zbar)';
end

% measurement update
Kk = Pxzbar*inv(Pzzbar);
x_hat = xbar + Kk*(z - zbar);
Pxx_hat = Pxxbar - Kk*Pzzbar*Kk';
