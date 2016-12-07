% exercise 2.20

x0 = [1; 1; 1; 0.1];

x = x0;
% use newton iteration to solve
for i=1:100
    [f, ja] = ex_2_20_funct(x);
    dx = ja\f;
    x = x-dx;
end

x
x(1) + x(2) + x(3)
x(4)
