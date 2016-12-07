% HW9 - Problem 2
% ASE 396
%
% D. Pederson


% read data file
phat = [10; pi/2];
Ppp = [0.02^2, 0; 0, (15*pi/180)^2];


                         
%% Part (a)
% compute mean and covariance via monte carlo
nmc = 100000;
p = repmat(phat', nmc,1) + randn(nmc,2)*sqrt(Ppp);
x = [p(:,1).*cos(p(:,2)), p(:,1).*sin(p(:,2))];
xbar_mc = [mean(x(:,1)), mean(x(:,2))]'
Pxx_mc = cov(x)

%% Part (b)
% mean and covariance via unscented transform
[xhat, Pxx] = UnscentedTransform(phat, Ppp, @(x) [x(1)*cos(x(2)); x(1)*sin(x(2))], 1, 3, 3-length(phat));

xhat
Pxx

%% Part (c)
% mean and covariance via taylor approx.
xbar_1st = [phat(1)*cos(phat(2)); phat(1)*sin(phat(2))]
dhdp = [cos(phat(2)), -phat(1)*sin(phat(2));
        sin(phat(2)), phat(1)*cos(phat(2))];
Pxx_1st = dhdp*Ppp*(dhdp')

%% Part (d)
% plot prob. 3sigma ellipsoid
fh = figure();
hold on
plot(x(:,1),x(:,2), 'k.')
xlabel('X')
ylabel('Y')
draw_cov_ellipse(gca, xbar_mc, Pxx_mc, 3)
draw_cov_ellipse(gca, xhat, Pxx, 3)
draw_cov_ellipse(gca, xbar_1st, Pxx_1st, 3)
legend('MC data points','MC','UT','1st Order')