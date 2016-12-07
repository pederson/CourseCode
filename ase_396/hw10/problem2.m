% HW10 - Problem 2
% ASE 396
%
% D. Pederson


% read data file
DAT = csvread('homework9_problem3.dat');
tdat = DAT(:,1);
rdat = DAT(:,2);
rrdat = DAT(:,3);
xtrudat = DAT(:,4);
xdtrudat = DAT(:,5);
z = [rdat, rrdat];
t = tdat;

% prior
x0 = [4.5;0.15];
%x0 = [4.9;0.25];
Pxx0 = [1000, 0; 0, 100];

% measurement error
R = [0.25^2, 0; 0, 0.1^2];


                         
%% Part (a)
x = x0;
Pxx = Pxx0;
Lyy = inv(Pxx0);
dy = zeros(size(x));
tprev = 0;
% setup Kalman filter
xhat = zeros(size(t));
xdhat = zeros(size(t));
resid = zeros(length(t),2);
Pk_hold = zeros(length(t), 2);
for k=1:length(t)
    
    [x, Lyy, dy] = InformationNonlinFilter(t(k)-tprev, z(k,:)', R, @(x)p2_measurement(x), ...
                                     @(t,x,u,v)p2_dynamics(t,x,u,v), ...
                                     @(t,x)zero_control(t,x), ... 
                                     @(t,x)p2_process_noise(t,x), ...
                                     x, Lyy, dy);
                                 
    
    dx = Lyy\dy;
    tprev = t(k);
    xhat(k) = x(1)+dx(1);
    xdhat(k) = x(2)+dx(2);
    Pxx = inv(Lyy);
    Pk_hold(k,1) = Pxx(1,1);
    Pk_hold(k,2) = Pxx(2,2);
    [hk, Hk] = p2_measurement(x);
    resid(k,:) = ((z(k,:)'-hk) - Hk*dx)';
end


% and the solution is:
x + dx
inv(Lyy)

figure()
hold on
plot(t, (xhat - xtrudat), 'k-')
%plot(t, xhat, 'b-')
%plot(t, xtrudat, 'k--')
plot(t, 3*sqrt(Pk_hold(:,1)), 'r--')
plot(t, -3*sqrt(Pk_hold(:,1)), 'r--')
title('With Process Noise')
xlabel('t')
ylabel('Position error')
%ylim([-2, 2])

% check the post-fit residuals
figure()
hold on
plot(t, resid(:,1), 'k-')
plot(t, resid(:,2), 'b-')
xlabel('t')
ylabel('residual')
title('Post-fit residual')
legend('\rho resid','\rho-dot resid')


