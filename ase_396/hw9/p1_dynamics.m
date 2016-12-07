function [xt, stm] = p1_dynamics(dt, x0)

omega = 2;

stm = [cos(omega*dt), 1/omega*sin(omega*dt);
       -omega*sin(omega*dt), cos(omega*dt)];
   
xt = stm*x0;