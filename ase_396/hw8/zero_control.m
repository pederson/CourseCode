function [vec, M] = zero_control(t, x0)

vec = zeros(size(x0));
M = zeros(length(x0), length(x0));

