% HW7 - Problem 1
% ASE 396
%
% D. Pederson

%% Part (b)
% read data file
DAT = csvread('homework7_problem1.dat');
tdat = DAT(:,1);
rdat = DAT(:,2);
z = rdat;
t = tdat;

% prior
x0 = [0.496; -1.527; 0.678; 0.734];
Pxx = 1e-6*eye(4);

% measurement error
R = 1e-7;

% setup GN
[xhat, Pxxhat] = GaussNewton(@(t,x)p1_dynamics(t,x), ... 
                             @(x)p1_measurement(x), ...
                             t, z, R, ...
                             x0, Pxx);

% and the solution is:
xhat
Pxxhat
   

% plot pre-fit and post-fit residuals
prefit = [];
for i=1:length(z)
    prefit = [prefit; (z(i,:)' - p1_measurement(x0))'];
end

postfit = [];
for i=1:length(z)
    postfit = [postfit; (z(i,:)' - p1_measurement(xhat))'];
end

figure()
hist(prefit(:,1))
title('Prefit - range')
%xlim([-3, 1])
figure()
hist(postfit(:,1))
title('Postfit - range')
%xlim([-3,1])