function [v, Q, gamma] = p1_process_noise(dt, x)

omega = 2;

gamma = [1/(omega^2)*(1-cos(omega*dt));
         1/omega*sin(omega*dt)];
   
v = 0;
Q = 1;