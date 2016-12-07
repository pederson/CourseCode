% Semester Project
% ASE 396
%
% D. Pederson


% read data file
DAT = csvread('ASE396_ProjectData.dat');
tdat = DAT(:,1);
rdat = DAT(:,2);
azdat = DAT(:,3);
eldat = DAT(:,4);
z = [rdat, azdat, eldat];
t = tdat;

% prior
x0 = [0; 0; 0; 0; 0; 0; 1e-3];
Pxx0 = diag([100, 100, 100, 10, 10, 10, 0.1]);

% measurement error
R = diag([0.5^2, (0.1*pi/180)^2, (0.2*pi/180)^2]);


                         
%% Part (2)
% running an EKF
x = x0;
Pxx = Pxx0;
tprev = 0;
state = zeros(length(t), length(x));
resid = zeros(length(t),length(z(1,:)));
innov = zeros(length(t),length(z(1,:)));
Pk_hold = zeros(length(t), length(x));
% setup Kalman filter
for k=1:length(t)
    
    [x, Pxx, inn] = ExtendedKalmanFilter(t(k)-tprev, z(k,:)', R, @(x,z)measurement(x,z), ...
                                     @(t,x,u,v)dynamics(t,x,u,v), ...
                                     @(t,x)zero_control(t,x), ... 
                                     @(t,x)process_noise(t,x, 1), ...
                                     x, Pxx);
                  
    tprev = t(k);
    state(k,:) = x';
    Pk_hold(k,:) = diag(Pxx)';
    [~,~,dzk] = measurement(x,z(k,:)');
    resid(k,:) = dzk';
    innov(k,:) = inn';
end


%% Part (3)
% tuning the filter
%a = logspace(1, 5, 5);
a = linspace(1, 100, 20);
%a = 100:50:100;

sigmarho = zeros(size(a));
sigmaalpha = zeros(size(a));
sigmabeta = zeros(size(a));
meanrho = zeros(size(a));
meanalpha = zeros(size(a));
meanbeta = zeros(size(a));
for j=1:length(a)
    x = x0;
    Pxx = Pxx0;
    tprev = 0;
    for k=1:length(t)

        [x, Pxx, inn] = ExtendedKalmanFilter(t(k)-tprev, z(k,:)', R, @(x,z)measurement(x,z), ...
                                         @(t,x,u,v)dynamics(t,x,u,v), ...
                                         @(t,x)zero_control(t,x), ... 
                                         @(t,x)process_noise(t,x,a(j)), ...
                                         x, Pxx);

        tprev = t(k);
        state(k,:) = x';
        Pk_hold(k,:) = diag(Pxx)';
        [~,~,dzk] = measurement(x,z(k,:)');
        resid(k,:) = dzk';
        innov(k,:) = inn';
    end
    sigmarho(j) = std(resid(:,1));
    sigmaalpha(j) = std(resid(:,2));
    sigmabeta(j) = std(resid(:,3));
    meanrho(j) = mean(resid(:,1));
    meanalpha(j) = mean(resid(:,2));
    meanbeta(j) = mean(resid(:,3));
end

figure()
hold on
plot(a, sigmarho, 'k-o')
plot(a, sigmaalpha, 'k-d')
plot(a, sigmabeta, 'k-s')
xlabel('a')
ylabel('Std. Dev.')
legend('\sigma_{\rho}','\sigma_{\alpha}','\sigma_{\beta}')

figure()
hold on
plot(a, abs(meanrho), 'k-o')
plot(a, abs(meanalpha), 'k-d')
plot(a, abs(meanbeta), 'k-s')
xlabel('a')
ylabel('Mean')
legend('\mu_{\rho}','\mu_{\alpha}','\mu_{\beta}')

%}

%% Part (4)
% plot the residuals
x = x0;
Pxx = Pxx0;
tprev = 0;
% setup Kalman filter
for k=1:length(t)
    
    [x, Pxx, inn] = ExtendedKalmanFilter(t(k)-tprev, z(k,:)', R, @(x,z)measurement(x,z), ...
                                     @(t,x,u,v)dynamics(t,x,u,v), ...
                                     @(t,x)zero_control(t,x), ... 
                                     @(t,x)process_noise(t,x, 20), ...
                                     x, Pxx);
                  
    tprev = t(k);
    state(k,:) = x';
    Pk_hold(k,:) = diag(Pxx)';
    [~,~,dzk] = measurement(x,z(k,:)');
    resid(k,:) = dzk';
    innov(k,:) = inn';
end


figure()
h=histogram(innov(:,1));
title('\rho Residuals')
xlabel('resid [m]')
ylabel('Count')
hold on
histogram(resid(:,1), h.BinEdges)
legend(['Innovation (\mu=',num2str(mean(innov(:,1))),'; \sigma=',num2str(std(innov(:,1))),')']...
      ,['Post-Fit (\mu=',num2str(mean(resid(:,1))),'; \sigma=',num2str(std(resid(:,1))),')'])


figure()
h=histogram(innov(:,2));
title('\alpha Residuals')
xlabel('resid [rad]')
ylabel('Count')
hold on
histogram(resid(:,2), h.BinEdges)
legend(['Innovation (\mu=',num2str(mean(innov(:,2))),'; \sigma=',num2str(std(innov(:,2))),')']...
      ,['Post-Fit (\mu=',num2str(mean(resid(:,2))),'; \sigma=',num2str(std(resid(:,2))),')'])


figure()
h=histogram(innov(:,3));
title('\beta Residuals')
xlabel('resid [rad]')
ylabel('Count')
hold on
histogram(resid(:,3), h.BinEdges)
legend(['Innovation (\mu=',num2str(mean(innov(:,3))),'; \sigma=',num2str(std(innov(:,3))),')']...
      ,['Post-Fit (\mu=',num2str(mean(resid(:,3))),'; \sigma=',num2str(std(resid(:,3))),')'])


% plot the estimated trajectory
figure()
scatter3(state(:,1), state(:,2), state(:,3),3,state(:,3))
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

% plot estimated individual states over time
figure()
subplot(3,1,1)
hold on
plot(t, state(:,1))
plot(t, state(:,1)+3*sqrt(Pk_hold(:,1)), 'r--')
plot(t, state(:,1)-3*sqrt(Pk_hold(:,1)), 'r--')
xlabel('T')
ylabel('x [m]')

subplot(3,1,2)
hold on
plot(t, state(:,2))
plot(t, state(:,2)+3*sqrt(Pk_hold(:,2)), 'r--')
plot(t, state(:,2)-3*sqrt(Pk_hold(:,2)), 'r--')
xlabel('T')
ylabel('y [m]')

subplot(3,1,3)
hold on
plot(t, state(:,3))
plot(t, state(:,3)+3*sqrt(Pk_hold(:,3)), 'r--')
plot(t, state(:,3)-3*sqrt(Pk_hold(:,3)), 'r--')
xlabel('T')
ylabel('z [m]')

figure()
subplot(3,1,1)
hold on
plot(t, state(:,4))
plot(t, state(:,4)+3*sqrt(Pk_hold(:,4)), 'r--')
plot(t, state(:,4)-3*sqrt(Pk_hold(:,4)), 'r--')
xlabel('T')
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')

subplot(3,1,2)
hold on
plot(t, state(:,5))
plot(t, state(:,5)+3*sqrt(Pk_hold(:,5)), 'r--')
plot(t, state(:,5)-3*sqrt(Pk_hold(:,5)), 'r--')
xlabel('T')
ylabel('$\dot{y}$ [m/s]','Interpreter','latex')

subplot(3,1,3)
hold on
plot(t, state(:,6))
plot(t, state(:,6)+3*sqrt(Pk_hold(:,6)), 'r--')
plot(t, state(:,6)-3*sqrt(Pk_hold(:,6)), 'r--')
xlabel('T')
ylabel('$\dot{z}$ [m/s]','Interpreter','latex')

figure()
hold on
plot(t, state(:,7))
plot(t, state(:,7)+3*sqrt(Pk_hold(:,7)), 'r--')
plot(t, state(:,7)-3*sqrt(Pk_hold(:,7)), 'r--')
xlabel('T')
ylabel('\omega [rad/s]')


% plot the variances
figure()
hold on
plot(t, Pk_hold(:,1))
plot(t, Pk_hold(:,2))
plot(t, Pk_hold(:,3))
plot(t, Pk_hold(:,4))
plot(t, Pk_hold(:,5))
plot(t, Pk_hold(:,6))
plot(t, Pk_hold(:,7))
xlabel('T')
ylabel('Covariance')
legend('P_{xx}','P_{yy}','P_{zz}',...
       'P_{$\dot{x}\dot{x}$}','P_{$\dot{y}\dot{y}$}','P_{$\dot{z}\dot{z}$}',...
       'P_{\omega\omega}')


% plot the individual measurements
figure()
subplot(3,1,1)
plot(t, z(:,1))
ylabel('\rho [m]')
xlabel('T')

subplot(3,1,2)
plot(t, z(:,2))
ylabel('\alpha [rad]')
xlabel('T')


subplot(3,1,3)
plot(t, z(:,3))
ylabel('\beta [rad]')
xlabel('T')



%% Part (5)
% get an a-priori using MLE of the first L data points
L = 50;
for s=1:L
    
end
