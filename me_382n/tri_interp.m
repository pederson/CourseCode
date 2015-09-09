function f = tri_interp(vertices, values, tri_points)
% Dylan Pederson
% This function creates an interpolant for a function given
% its values at the vertices of the triangle
numpts = length(tri_points)/2;

% coefficients for the interpolant
gam = [values(1); values(2)-values(1); values(3)-values(1)];

% coefficients for the x coordinate
a = [vertices(1); vertices(3)-vertices(1); vertices(5)-vertices(1)];

% coefficients for the y coordinate
b = [vertices(2); vertices(4)-vertices(2); vertices(6)-vertices(2)];

% ut and vt and fi
f = zeros(numpts, 1);
lhm = [a(2), a(3);b(2), b(3)];
invlhm = inv(lhm);
for i=1:numpts
    rhs = [tri_points(2*i-1)-a(1);tri_points(2*i)-b(1)];
    uv = invlhm*rhs;
    f(i) = gam(1) + gam(2)*uv(1) + gam(3)*uv(2);
end

end