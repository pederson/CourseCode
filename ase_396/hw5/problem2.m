% HW 5 - problem 2
%
% ASE 396
%
% D. Pederson


% read data from file
D=dlmread('homework5_problem2.dat');
t = D(:,1);
xs = D(:,2);
ys = D(:,3);
zi = D(:,4);

% plot prefit residuals
eps_pre = zi - sqrt(xs.^2 + ys.^2);
figure()
hist(eps_pre)
title('Prefit residuals')

% use nonlinear least squares (GN) 
R = 2*eye(length(zi));
Rinv = 0.5*eye(length(zi));
xstar = zeros(2,1);
Pxx = 100*eye(2);
dxbar = zeros(2,1);
dz = eps_pre;
for n=1:3
    % calculate new H
  H = [xstar(1)-xs, xstar(2)-ys];  
    for i=1:length(xs)
        H(i) = H(i)/sqrt((xstar(1)-xs(i))^2 + (xstar(2)-ys(i))^2);
    end
    
    % solve least squares problem
    dxhat = inv(H'*Rinv*H + inv(Pxx))*(H'*Rinv*dz + inv(Pxx)*dxbar);
    
    % update 
    dxbar = dxbar - dxhat;
    xstar = xstar + dxhat;
end

% calculate the post-fit residuals
eps_post = zi - sqrt((xstar(1)-xs).^2 + (xstar(2)-ys).^2);
figure()
hist(eps_post)
title('Postfit residuals')