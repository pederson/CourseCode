function xdot = central_force(x)

xdot = [x(3); x(4); -x(1)./(x(1).^2 +x(2).^2)^(3/2); -x(2)./(x(1).^2 +x(2).^2)^(3/2)];
