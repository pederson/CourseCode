function [] = hw5_1(A,b, epsilon, xstar)

%% part A
% compute the reduced SVD and plot singular values
[U, S, V] = svd(A, 0);
sing_vals = diag(A);
semilogy(sing_vals);

%% part B
% explain whether A is rank deficient (or nearly)


%% part C
% test svd, qr w/ pivot, qr w/o pivot, and rangeA.m function
% to find the best rank-r approximation to Range(A)

% svd

% qr w/ pivot

% qr w/o pivot

% rangeA.m


%% part D
% use truncated SVD to solve least squares

% (i)
% plot relative cond, pertub A

% (ii)
% plot relative cond, pertub b

% (iii)
% plot the norm |xr|/|xstart| as fn(r)

% (iv)
% plot residual as fn(r)

% (v)
% plot the error as fn(r)

% (vi)
% plot the residual v error for each r using loglog


%% part E
% repeat part D but with regularized svd (beta = sigma_r)

