% HW9 - Problem 3
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
%x0 = [4.5;0.15];
x0 = [4.9;0.25];
Pxx0 = [1000, 0; 0, 100];

% measurement error
R = [0.25^2, 0; 0, 0.1^2];


                         
%% Part (a)
x = x0;
Pxx = Pxx0;
tprev = 0;
% setup Kalman filter
xhat = zeros(size(t));
xdhat = zeros(size(t));
resid = zeros(length(t),2);
Pk_hold = zeros(length(t), 2);
for k=1:length(t)
    
    [x, Pxx] = UnscentedKalmanFilter(t(k)-tprev, z(k,:)', R, @(x)p3_measurement(x), ...
                                     @(t,x,u,v)p3_dynamics(t,x,u,v), ...
                                     @(t,x)zero_control(t,x), ... 
                                     @(t,x)p3_process_noise(t,x), ...
                                     x, Pxx);
                                 
    
    tprev = t(k);
    xhat(k) = x(1);
    xdhat(k) = x(2);
    Pk_hold(k,1) = Pxx(1,1);
    Pk_hold(k,2) = Pxx(2,2);
    resid(k,:) = (z(k,:)'-p3_measurement(x))';
end


% and the solution is:
x
Pxx

figure()
hold on
plot(t, (xhat - xtrudat), 'k-')
% plot(t, xtrudat, 'k--')
plot(t, 3*sqrt(Pk_hold(:,1)), 'r--')
plot(t, -3*sqrt(Pk_hold(:,1)), 'r--')
title('With Process Noise')
xlabel('t')
ylabel('Position error')
ylim([-2, 2])

% check the post-fit residuals
figure()
hold on
plot(t, resid(:,1), 'k-')
plot(t, resid(:,2), 'b-')
xlabel('t')
ylabel('residual')
title('Post-fit residual')
legend('\rho resid','\rho-dot resid')


%{
Estimated state stays within the bounds of the estimated covariance... so
yes it is consistent with the truth
%}
