% HW 5 - problem 1
%
% ASE 396
%
% D. Pederson

z = [0, 0.75, -0.8]';
H = [1, 1; 0, 1; 1, 0];
xbar = zeros(2,1);
Pxx = eye(2);
R = 0.5*eye(3);
Rinv = 2*eye(3);


%% batch estimators

% MAP estimate
K = Pxx*H'*inv(H*Pxx*H'+R);
x_map = xbar + K*(z-H*xbar)
Pxxz_map = Pxx - K*H*Pxx


% LS with a priori
Pxxz_ls = inv(H'*Rinv*H + inv(Pxx))
x_ls = Pxxz_ls*(H'*Rinv*z+inv(Pxx)*xbar)



%% sequential estimation part
R = 0.5;
Rinv = 2;

% first
z = 0;
H = [1,1];
xbar = zeros(2,1);
Pxx = eye(2);
K = Pxx*H'*inv(H*Pxx*H'+R);
xbar = xbar + K*(z-H*xbar);
Pxx = Pxx - K*H*Pxx;

% second
z = 0.75;
H = [0,1];
K = Pxx*H'*inv(H*Pxx*H'+R);
xbar = xbar + K*(z-H*xbar);
Pxx = Pxx - K*H*Pxx;

% third
z = -0.8;
H = [1,0];
K = Pxx*H'*inv(H*Pxx*H'+R);
xbar = xbar + K*(z-H*xbar)
Pxx = Pxx - K*H*Pxx

