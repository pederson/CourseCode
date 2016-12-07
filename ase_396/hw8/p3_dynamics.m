function [xdot] = p3_dynamics(x)


stm = x(5:end);
stm = reshape(stm, [4,4]);

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     -2*x(2)/(x(1)+x(2))^3, 1/(x(1)+x(2))^2 - 2*x(2)/(x(1)+x(2))^3, 0, 0;
     1/(x(1)+x(2))^2 - 2*x(1)/(x(1)+x(2))^3, -2*x(1)/(x(1)+x(2))^3, 0, 0];
 
xdot = [x(3); x(4); x(2)/(x(1)+x(2))^2; x(1)/(x(1)+x(2))^2];
stmdot = A*stm;


xdot = [xdot; reshape(stmdot, [16,1])];