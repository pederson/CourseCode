function [xt, stm] = p2_dynamics(t, x0)

omega = 2;

stm = [cos(omega*t), 1/omega*sin(omega*t);
       -omega*sin(omega*t), cos(omega*t)];
   
xt = stm*x0;