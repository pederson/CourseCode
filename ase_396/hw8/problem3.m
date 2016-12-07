% HW8 - Problem 3
% ASE 396
%
% D. Pederson

% prior
x0 = [0.3; -1.3; -0.433; 0.342];

% set STM IC
STM_0 = eye(4);
STM_0 = reshape(STM_0, [16, 1]);

x0 = [x0; STM_0];

                 
%% Part (a)
ti = 100*(0:10);

% numerically integrate
opts = odeset('RelTol',1e-15,'AbsTol',ones(20,1)*1e-15);
[~, xnum] = ode45(@(t,x) p3_dynamics(x), ti, x0, opts);

% here is x(t1)
x1a = xnum(2,1:4)

% here is x(t10)
x10a = xnum(11,1:4)

%% Part (b)

% perturb the initial data
x0 = [0.3; -1.3; -0.433; 0.342];
dx0 = 1e-4*[1; 1; -1; -1];
x0 = x0 + dx0;
x0 = [x0; STM_0];

% time
ti = 100*(0:10);

% numerically integrate
opts = odeset('RelTol',1e-15,'AbsTol',ones(20,1)*1e-15);
[~, xnum] = ode45(@(t,x) p3_dynamics(x), ti, x0, opts);

% here is x(t1)
x1b = xnum(2,1:4)

% here is STM(t1)
STM1b = reshape(xnum(2,5:end),[4,4])

% here is x(t10)
x10b = xnum(11,1:4)

% here is STM(t10)
STM10b = reshape(xnum(11,5:end), [4,4])

%% Part (c)

% at time t1
dx1_1 = x1b - x1a
dx1_2 = (STM1b*dx0)'
delta1 = dx1_1 - dx1_2

% at time t10
dx10_1 = x10b - x10a
dx10_2 = (STM10b*dx0)'
delta10 = dx10_1 - dx10_2


%{
 The two are fairly close... this differences are on the order of 10^-5
 . This indicates that the propagated STM is pretty accurate, but the
 accuracy may decrease over longer time spans
%}