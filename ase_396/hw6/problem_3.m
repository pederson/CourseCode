% HW6 - Problem 3 
% ASE 396
%
% D. Pederson

tvec = 0:10;

%% Part (a)
% initial condition
x_0 = [1; 0; 0; 1];

% the analytical solution
% x_true = A*cos(sqrt(k_m)*tvec+phi);

% integrate
opts = odeset('RelTol',1e-7,'AbsTol',1e-5);
[~, x_num] = ode45(@(t,x) central_force(x), tvec, x_0, opts);

% plot the result
figure()
hold on
% plot(tvec, x_true, 'k-')
% plot(tvec, x_num(:,1), 'kd')
plot(x_num(:,1), x_num(:,2), 'kd')
xlabel('X (DU)')
ylabel('Y (DU)')
% legend('Analytical','Numerical')
title(['Central Force - Noise Free'])

% here's x(t1)
x1 = x_num(2,:)

% here's x(t10)
x10 = x_num(11,:)

%% Part (b)

% set new IC
dx_0 = 1e-6*[-1; 1; -1; -1];
x_0 = x_0 + dx_0;

% set STM IC
STM_0 = eye(4);
STM_0 = reshape(STM_0, [16, 1]);

x_0 = [x_0; STM_0];

% integrate simultaneously
opts = odeset('RelTol',1e-7,'AbsTol',1e-5);
[~, x_num_b] = ode45(@(t,x) central_force_STM(x), tvec, x_0, opts);

% plot the result
figure()
hold on
% plot(tvec, x_true, 'k-')
% plot(tvec, x_num(:,1), 'kd')
plot(x_num(:,1), x_num(:,2), 'kd')
plot(x_num_b(:,1), x_num_b(:,2), 'rd')
xlabel('X (DU)')
ylabel('Y (DU)')
% legend('Analytical','Numerical')
title(['Central Force - Noise Free'])

% here's x(t1)
xstar1 = x_num_b(2,1:4)

% here's x(t10)
xstar10 = x_num_b(11,1:4)

% here's STM(t1)
STM1 = reshape(x_num_b(2, 5:end), [4,4])

% here's STM(t10)
STM10 = reshape(x_num_b(11, 5:end), [4,4])



%% Part (c)

% report deltas at t1
dx1_1 = xstar1 - x1
dx2_1 = STM1*dx_0


% report deltas at t10
dx1_10 = xstar10 - x10
dx2_10 = STM10*dx_0


%{
 There's about an order of magnitude difference
 in the two delta vectors at t10. The delta as propagated 
 the STM is an order of magnitude SMALLER than the straight integrated
 x values. This tells me that using the STM to propagate is actually
 more accurate, and less sensitive to noise/perturbations.
 

%}