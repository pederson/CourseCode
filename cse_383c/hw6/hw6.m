% hw6
%
% D. Pederson
%
% Driver program the does the following:
% (i) calls tridiag
% (ii) calls qralg to get one eigenvalue
% (iii) calls qralg with smaller matrix to get another eigenvalue
% (iv) continue for all eigenvalues


A = hilb(4);
%A = diag(15:-1:1) + ones(15, 15);

% call tridiag
T = tridiag(A);

% get one eigenvalue with qralg
