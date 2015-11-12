function [] = hw5_2()

% load the matrix
o = matrixfree;
m = o.m;
n = o.n;
fullrank = min([m,n]);

disp(['[m,n] = ',num2str(m), ',', num2str(n)])

allranks = 1:1:100;%200;
resid = zeros(size(allranks));

for i = 1:length(allranks)
    targetrank = allranks(i);
    
    % form orthonormal basis
    Z = rangeA(@(x)o.A(x), @(x)o.At(x), m, n, targetrank);

    % form B = Z'A
    B = zeros(targetrank, n);
    size(B)
    for j=1:targetrank
        zj = Z(:,j);
        bj = o.At(zj);
        bj = bj';
        B(j, :) = bj;
    end

    % form the QR factorization of B
    [W,T] = qr(B);

    % now the decomp. is Z T' W'
    % the decomp. is Z W T
    % solve least squares system
    y = Z'*o.b;
    q = W'*y;
    x = T\q;

    % do some error analysis
    % residual
    resid(i) = norm(o.b - o.A(x))/norm(o.b);
end

plot(allranks, resid, 'k')

