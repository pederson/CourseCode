function [h, H] = p1_measurement(x)

height = 5.4;
rho = sqrt(x(1)^2 + height^2);

h = [rho;
     x(1)*x(2)/rho];
 
H = [x(1)/rho                        ,    0;
     x(2)/rho - x(1)^2*x(2)/rho^3    ,    x(1)/rho];