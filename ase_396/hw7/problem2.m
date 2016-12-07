% HW7 - Problem 2
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
x0 = [5;1];
Pxx = [1000, 0; 0, 100];

% measurement error
R = [0.25^2, 0; 0, 0.1^2];

%% Part (c)
% setup GN
[xhat, Pxxhat] = GaussNewton(@(t,x)p2_dynamics(t,x), ... 
                             @(x)p2_measurement(x), ...
                             t, z, R, ...
                             x0, Pxx);
 
                         
% and the solution is:
xhat
Pxxhat
       

%% Part (d)
% plot pre-fit and post-fit residuals
prefit = [];
for i=1:length(z)
    prefit = [prefit; (z(i,:)' - p2_measurement(x0))'];
end

postfit = [];
for i=1:length(z)
    postfit = [postfit; (z(i,:)' - p2_measurement(xhat))'];
end

figure()
hist(prefit(:,1))
title('Prefit - range')
xlim([-3,1])
figure()
hist(postfit(:,1))
title('Postfit - range')
xlim([-3,1])

figure()
hist(prefit(:,2))
title('Prefit - rangerate')
xlim([-4,3])
figure()
hist(postfit(:,2))
title('Postfit - rangerate')
xlim([-4, 3])

