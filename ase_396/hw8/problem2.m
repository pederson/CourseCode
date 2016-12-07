% HW8 - Problem 2
% ASE 396
%
% D. Pederson


% read data file
DAT = csvread('homework8_problem2.dat');
tdat = DAT(:,1);
rdat = DAT(:,2);
rrdat = DAT(:,3);
xtrudat = DAT(:,4);
xdtrudat = DAT(:,5);
z = [rdat, rrdat];
t = tdat;

% prior
x0 = [4.5;0.15];
Pxx0 = [1000, 0; 0, 100];

% measurement error
R = [0.25^2, 0; 0, 0.1^2];


                         
%% Part (a)
x = x0;
dx = 0*x;
Pxx = Pxx0;
tprev = 0;
% setup Kalman filter
xhat = zeros(size(t));
xdhat = zeros(size(t));
Pk_hold = zeros(length(t), 2);
for k=1:length(t)
    
    [x, Pxx, dx] = ClassicalKalmanNonlinFilter(t(k)-tprev, z(k,:)', R, @(x)p1_measurement(x), ...
                                     @(t,x)p1_dynamics(t,x), ...
                                     @(t,x)zero_control(t,x), ... 
                                     @(t,x)zero_noise(t,x), ...
                                     x, Pxx, dx);
    tprev = t(k);
    xhat(k) = x(1);
    xdhat(k) = x(2);
    Pk_hold(k,1) = Pxx(1,1);
    Pk_hold(k,2) = Pxx(2,2);
end


% and the solution is:
x
Pxx
dx

figure()
hold on
plot(t, (xhat - xtrudat), 'k--')
% plot(t, xtrudat, 'k--')
plot(t, 3*sqrt(Pk_hold(:,1)), 'r--')
plot(t, -3*sqrt(Pk_hold(:,1)), 'r--')
title('Without Process Noise')
xlabel('t')
ylabel('Position error')
ylim([-2, 2])

%{
figure()
hold on
plot(t, xdhat -xdtrudat, 'k--')
% plot(t, xdtrudat, 'k--')
plot(t, 3*sqrt(Pk_hold(:,2)), 'r--')
%}


%% Part (b)

%{ 
See attachment for handwritten derivation
%}

%% Part (c)
x = x0;
dx = 0*x;
Pxx = Pxx0;
tprev = 0;
% setup Kalman filter
for k=1:length(t)
    
    [x, Pxx] = ClassicalKalmanNonlinFilter(t(k)-tprev, z(k,:)', R, @(x)p1_measurement(x), ...
                                     @(t,x)p1_dynamics(t,x), ...
                                     @(t,x)zero_control(t,x), ... 
                                     @(t,x)p1_process_noise(t,x), ...
                                     x, Pxx, dx);
                                 
    tprev = t(k);
    xhat(k) = x(1);
    xdhat(k) = x(2);
    Pk_hold(k,1) = Pxx(1,1);
    Pk_hold(k,2) = Pxx(2,2);
end


% and the solution is:
x
Pxx


figure()
hold on
plot(t, (xhat - xtrudat), 'k--')
% plot(t, xtrudat, 'k--')
plot(t, 3*sqrt(Pk_hold(:,1)), 'r--')
plot(t, -3*sqrt(Pk_hold(:,1)), 'r--')
title('With Process Noise')
xlabel('t')
ylabel('Position error')


%{
Estimated state stays within the bounds of the estimated covariance... so
yes it is consistent with the truth
%}
