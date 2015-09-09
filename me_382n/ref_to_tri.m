function [tri_points] = ref_to_tri(vertices, ref_points)
% Dylan Pederson
% This function maps points in the reference
% triangle to points on the given triangle
num_pts = length(ref_points)/2;

% coefficients for the x coordinate
a = [vertices(1); vertices(3)-vertices(1); vertices(5)-vertices(1)];

% coefficients for the y coordinate
b = [vertices(2); vertices(4)-vertices(2); vertices(6)-vertices(2)];

% map the reference points back to normal x, y points
hold_mat = reshape(ref_points, num_pts, 2)';
ref_mat = [ones(num_pts, 1), hold_mat];
x_points = ref_mat*a;
y_points = ref_mat*b;

%{
figure()
verts = reshape(vertices, 2, 3)';
hold on
plot(x_points, y_points, 'b*')
plot(hold_mat(:,1), hold_mat(:,2), 'k.')
plot([verts(:,1); verts(1,1)], [verts(:,2); verts(1,2)], 'b')
plot([0, 1, 0, 0], [0, 0, 1, 0], 'k')
%}

tri_points = reshape([x_points, y_points]', num_pts*2, 1);
end