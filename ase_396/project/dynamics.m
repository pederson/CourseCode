function [xt, stm] = dynamics(dt, x, u, v)

if dt == 0
    xt = x;
    stm = eye(length(x));
    return;
end

w = x(end);

dxdw = -x(5)*dt*cos(dt*w) - x(4)*dt*sin(dt*w);
dydw = x(4)*dt*cos(dt*w) - x(5)*dt*sin(dt*w);

if w == 0
    stm = [1, 0, 0, dt, 0, 0, 0;
           0, 1, 0, 0, dt, 0, 0;
           0, 0, 1, 0, 0, dt, 0;
           0, 0, 0, 1, 0, 0, dxdw;
           0, 0, 0, 0, 1, 0, dydw;
           0, 0, 0, 0, 0, 1, 0;
           0, 0, 0, 0, 0, 0, 1];
       
    gamma = [dt^2/2, 0, 0, 0;
             0, dt^2/2, 0, 0;
             0, 0, dt^2/2, 0;
             dt, 0, 0, 0;
             0, dt, 0, 0;
             0, 0, dt, 0;
             0, 0, 0, dt];
else
    stm = [1, 0, 0, sin(w*dt)/w, -(1-cos(w*dt))/w, 0, 0;
           0, 1, 0, (1-cos(w*dt))/w, sin(w*dt)/w, 0, 0;
           0, 0, 1, 0, 0, dt, 0;
           0, 0, 0, cos(w*dt), -sin(w*dt), 0, dxdw;
           0, 0, 0, sin(w*dt), cos(w*dt), 0, dydw;
           0, 0, 0, 0, 0, 1, 0;
           0, 0, 0, 0, 0, 0, 1];
       
       
    gamma = [(1-cos(w*dt))/w^2, -(w*dt-sin(w*dt))/w^2, 0, 0;
             (w*dt-sin(w*dt))/w^2, (1-cos(w*dt))/w^2, 0, 0;
             0, 0, dt^2/2, 0;
             sin(w*dt)/w, -(1-cos(w*dt))/w, 0, 0;
             (1-cos(w*dt))/w, sin(w*dt)/w, 0, 0;
             0, 0, dt, 0;
             0, 0, 0, dt];

end

xt = stm*x + gamma*v;
end