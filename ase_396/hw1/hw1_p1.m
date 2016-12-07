%% HW1 - Problem 1
% Differential equation:
%       d/dt(dx/dt) = -k/m x
%
% turning the problem into a matrix system
%
%       d/dt(x) = v
%       d/dt(v) = -k/m x


k_m = 1;
phi = pi/3;
A = 1.34;
tvec = 0:20;

% initial condition
x_0 = [A*cos(phi); -A*sqrt(k_m)*sin(phi)];

% the analytical solution
x_true = A*cos(sqrt(k_m)*tvec+phi);

% form the B matrix
B = [0 1; -k_m 0];

% integrate
opts = odeset('RelTol',1e-7,'AbsTol',1e-5);
[~, x_num] = ode45(@(t,x) harmonic_oscillator(B, x), tvec, x_0, opts);

% plot the result
figure()
hold on
plot(tvec, x_true, 'k-')
plot(tvec, x_num(:,1), 'kd')
xlabel('Time (TU)')
ylabel('X (DU)')
legend('Analytical','Numerical')
title(['Harmonic Oscillator',' A=',num2str(A),' k/m=',num2str(k_m),' \phi=',num2str(phi)])

% plot the numerical error over time
err = x_true-x_num(:,1)';
figure()
plot(tvec, err,'k')
xlabel('Time (TU)')
ylabel('Error (DU)')
title('Numerical Error')