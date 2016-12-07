function [h, H, dz] = measurement(x, z)


xs = 1000;
ys = 1000;
zs = 20;

xrel = x(1) - xs;
yrel = x(2) - ys;
zrel = x(3) - zs;
rho = sqrt(xrel^2 + yrel^2 + zrel^2);

% the measurement function
h = [rho; atan2(xrel, yrel); asin(zrel/rho)];

% measurement jacobian
H = [xrel/rho, yrel/rho, zrel/rho, 0, 0, 0, 0;
     yrel/(xrel^2 + yrel^2), -xrel/(xrel^2 + yrel^2), 0, 0, 0, 0, 0;
     -xrel*zrel/(rho^2*sqrt(rho^2 - zrel^2)), ...
     -yrel*zrel/(rho^2*sqrt(rho^2 - zrel^2)), ...
     (1-zrel^2/rho^2)/sqrt(rho^2 - zrel^2), 0, 0, 0, 0];

% also output the residual
zp = z;
if z(2) < 0
    zp(2) = zp(2) + 2*pi;
end
hp = h;
if h(2) < 0
    hp(2) = hp(2) + 2*pi;
end
dz = zp - hp;