function [xt, stm] = p1_dynamics(t, x0)

if t==0
    xt = x0;
    stm = eye(4);
    return;
end

dt = 1;
tvec = 0:dt:t;

STM_0 = eye(4);
STM_0 = reshape(STM_0, [16, 1]);

x_0 = [x0; STM_0];

% integrate simultaneously
opts = odeset('RelTol',1e-7,'AbsTol',1e-5);
[~, x_num] = ode45(@(t,x) central_force_STM(x), tvec, x_0, opts);

% here's x(t)
xt = x_num(end,1:4);

% here's STM(t)
stm = reshape(x_num(end, 5:end), [4,4]);