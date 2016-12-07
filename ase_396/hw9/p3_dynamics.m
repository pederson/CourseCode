function [xt] = p3_dynamics(dt, x0, u, v)

if dt == 0
    xt = x0;
    return;
end

%{
gamma = [1/(omega^2)*(1-cos(omega*dt));
         1/omega*sin(omega*dt)];

stm = [cos(omega*dt), 1/omega*sin(omega*dt);
       -omega*sin(omega*dt), cos(omega*dt)];
%}


tvec = 0:dt/10:dt;
opts = odeset('RelTol',1e-13,'AbsTol',ones(length(x0),1)*1e-15);
[~, x_num] = ode45(@(t,x) forcing(t, x, u, v), tvec, x0, opts);

xt = x_num(end,:)';
end



function xdot = forcing(t, x, u, v)

    omega = 2;
    xdot = [x(2); -omega^2*x(1)+v];

end