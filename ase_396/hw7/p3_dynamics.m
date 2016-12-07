function [xt, stm] = p3_dynamics(t, x0)

stm = [1, 1; 0, 1];
xt = stm*x0;