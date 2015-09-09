function interp_points = quad_interp(points, values, eval_points)
% Dylan Pederson

% quadratic interpolator that assumes the input points are 
% in ascending order and requires the number of input points
% to be odd

% make sure that points is the same size as values
if length(points) ~= length(values)
    error('points and values must be the same length')
end

% make sure evaluation points are within the bounds of the input points
if nnz(eval_points > max(points)) > 0 || nnz(eval_points < min(points)) > 0
    error('some evaluation points are out of bounds...cannot extrapolate');
end

% make sure there is an odd number of input points
if mod(length(points),2) == 0
    error('must input an odd number of points')
end

interp_points = zeros(size(eval_points));
i2 = zeros(size(eval_points));
% for each evaluation point
valid = true(size(eval_points));
for i=1:length(eval_points)

    % find the correct 3 points to interpolate
    pos_finder = points - eval_points(i);
    
    % check if the point is exactly one of the inputs
    if nnz(pos_finder == 0) ~= 0
        i2(i) = find(pos_finder == 0);
        valid(i) = false;
        continue;
    end

    % find where it goes from positive to negative
    position = find(pos_finder > 0);
    position = position(1);
    
    % find the middle index
    if mod(position,2) == 0
        i2(i) = position;
    else
        i2(i) = position - 1;
    end
end

% calculate the other two indices
i1 = i2-1;
i3 = i2+1;

% calculate interpolation points
interp_points(valid) = interp_quad_three(points, values,...
                                        i1(valid), i2(valid), i3(valid), eval_points(valid));
interp_points(~valid) = values(i2(~valid));

end

% This function does the actual interpolation
% it is set up so that it can be called in a vectorized way
function interp_val = interp_quad_three(points, values, i1, i2, i3, eval_point)
    x1 = points(i1);
    x2 = points(i2);
    x3 = points(i3);
    
    f1 = values(i1);
    f2 = values(i2);
    f3 = values(i3);
    
    interp_val = f1.*(eval_point - x2).*(eval_point - x3)./((x1-x2).*(x1-x3)) + ...
                 f2.*(eval_point - x1).*(eval_point - x3)./((x2-x3).*(x2-x1)) + ...
                 f3.*(eval_point - x1).*(eval_point - x2)./((x3-x2).*(x3-x1));
             
end