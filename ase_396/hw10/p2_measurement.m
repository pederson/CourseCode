function [h, H] = p2_measurement(x)

height = 135;
x0 = -100;
rho = sqrt((x(1)-x0)^2 + height^2);

h = [rho;
     (x(1)-x0)*x(2)/rho];
 
H = [(x(1)-x0)/rho                        ,    0;
     x(2)/rho - (x(1)-x0)^2*x(2)/rho^3    ,    (x(1)-x0)/rho];