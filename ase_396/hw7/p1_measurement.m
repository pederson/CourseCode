function [h, H] = p1_measurement(x)

rho = sqrt(x(1)^2 + x(2)^2);

h = rho;
 
H = [x(1)/rho, x(2)/rho, 0, 0];