function [] = hw5_2()

o = matrixfree;
m = o.m;
n = o.n;
targetrank = 10;

% form orthonormal basis
Z = rangeA(@(x)o.A(x), @(x)o.At(x), m, n, targetrank);

% form B = Z'A
B = zeros(targetrank, n);
for i=1:targetrank
    zi = Z(:,i);
    bi = o.At(zi);
    bi = bi';
    B(i, :) = bi;
end

% form the QR factorization of B
[W,T] = qr(B');

% now the decomp. is Z T' W'
% solve least squares system
