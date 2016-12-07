function xdot = central_force_STM(x)
% here, the x vector contains the STM appended

rsq = x(1).^2 +x(2).^2;

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     -1/rsq^(3/2)+3*x(1)^2/rsq^(5/2), 3*x(1)*x(2)/rsq^(5/2), 0, 0;
     3*x(1)*x(2)/rsq^(5/2), -1/rsq^(3/2)+3*x(2)^2/rsq^(5/2), 0, 0];

xdot = [x(3); 
        x(4); 
        -x(1)./rsq^(3/2); 
        -x(2)./rsq^(3/2);];
    
STM = reshape(x(5:end), [4, 4]);

STM_dot = A*STM;

xdot = [xdot; reshape(STM_dot, [16, 1])];
