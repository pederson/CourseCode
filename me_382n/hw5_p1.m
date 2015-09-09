% hw5 problem 1

% Solving Ax=b

A = gallery('prolate',32);
b = ones(32,1);

% solve using LU
[L,U,P] = lu(A);
x_lu = U\(L\(P*b));

% solve using pcg
x_pcg = pcg(A, b);

% add noise and solve using LU
b_hat = b+0.001*rand(32,1);
x_lu_hat = U\(L\(P*b_hat));

% calculate the errors
err_pcg = norm(x_pcg-x_lu)/norm(x_lu)
err_hat = norm(x_lu_hat-x_lu)/norm(x_lu)

