function [vec, M, N] = zero_noise(t, x0)

vec = zeros(size(x0));
M = zeros(length(x0), length(x0));
N = M;