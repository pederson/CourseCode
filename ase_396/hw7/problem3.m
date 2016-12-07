% HW7 - Problem 3
% ASE 396
%
% D. Pederson


% read data file
z = [3, 2;
     1, 7];
t = [1; 2];

% prior
x0 = [2; 3];
Pxx0 = eye(2);

% measurement error
R = [1, 0; 0, 2];

%% Part (a)
x = x0;
Pxx = Pxx0;
% setup Kalman filter
for k=1:length(t)
    
    [x, Pxx] = ClassicalKalmanFilter(@(t,x)p3_dynamics(t,x), ... 
                             @(x)p3_measurement(x), ...
                             t, z(k,:)', R, ...
                             x, Pxx);
end
 
                         
% and the solution is:
x
Pxx


%% Part (b)


% estimate the best x0 using G-N
[xhat, Pxxhat] = GaussNewton(@(t,x)p3_dynamics(t,x), ... 
                             @(x)p3_measurement(x), ...
                             t, z, R, ...
                             x0, Pxx);

% and the solution is:
xhat
Pxxhat