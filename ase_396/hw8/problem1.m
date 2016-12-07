% HW8 - Problem 1
% ASE 396
%
% D. Pederson


% read data file
DAT = csvread('homework7_problem2.dat');
tdat = DAT(:,1);
rdat = DAT(:,2);
rrdat = DAT(:,3);
z = [rdat, rrdat];
t = tdat;

% prior
x0 = [4.5;0.15];
Pxx0 = [1000, 0; 0, 100];

% measurement error
R = [0.25^2, 0; 0, 0.1^2];


                         
%% Part (a)
x = x0;
Pxx = Pxx0;
dx = 0*x;
tprev = 0;
% setup Kalman filter
Pk_hold = zeros(length(t), 2);
correl= zeros(length(t),1);
for k=1:length(t)
    
    [x, Pxx, dx] = ClassicalKalmanNonlinFilter(t(k)-tprev, z(k,:)', R, @(x)p1_measurement(x), ...
                                     @(t,x)p1_dynamics(t,x), ...
                                     @(t,x)zero_control(t,x), ... 
                                     @(t,x)zero_noise(t,x), ...
                                     x, Pxx, dx);
    tprev = t(k);
    Pk_hold(k,1) = Pxx(1,1);
    Pk_hold(k,2) = Pxx(2,2);
    correl(k) = Pxx(2,1)/(sqrt(Pxx(1,1))*sqrt(Pxx(2,2)));
end


% and the solution is:
x
Pxx


%% Part (b)

[xhat, Pxxhat] = GaussNewton(@(t,x)p1_dynamics(t,x), ... 
                             @(x)p1_measurement(x), ...
                             t, z, R, ...
                             x0, Pxx0, 1);

% and the solution is:
xhat
Pxxhat

dxhat_0 = xhat - x0;
dxhat_f = p1_dynamics(t(end), dxhat_0)
dxhat_KF = dx


%% Part (c)
% diagonal elements of Pk
semilogy(t, Pk_hold(:,1), 'r-.')
hold on
semilogy(t, Pk_hold(:,2), 'k--')
legend('P_{xx}','P_{vv}')
xlabel('t')

%{
Covariance decreases over time as expected... filter saturation when there
is zero process noise 
%}


%% Part (d)
% correlation coefficient over time
figure()
plot(t, correl, 'ko-')
xlabel('t')
ylabel('Correlation Coefficient')

%{
No issues here, as long as the correlation stays comfortably far away from
1... there is one instance where it is close to 1
%}
