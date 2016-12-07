function [v, Q, gamma] = process_noise(dt, x, a)

w = x(end);

if w == 0
    gamma = [dt^2/2, 0, 0, 0;
             0, dt^2/2, 0, 0;
             0, 0, dt^2/2, 0;
             dt, 0, 0, 0;
             0, dt, 0, 0;
             0, 0, dt, 0;
             0, 0, 0, dt];
else
    gamma = [(1-cos(w*dt))/w^2, -(w*dt-sin(w*dt))/w^2, 0, 0;
             (w*dt-sin(w*dt))/w^2, (1-cos(w*dt))/w^2, 0, 0;
             0, 0, dt^2/2, 0;
             sin(w*dt)/w, -(1-cos(w*dt))/w, 0, 0;
             (1-cos(w*dt))/w, sin(w*dt)/w, 0, 0;
             0, 0, dt, 0;
             0, 0, 0, dt];
end
   
v = zeros(4,1);
Q = a*diag([1e-3, 1e-3, 1e-3, 1e-7]);