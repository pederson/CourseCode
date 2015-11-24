% hw6
%
% D. Pederson
%
% Driver program that does the following:
% (i) calls tridiag
% (ii) calls qralg to get one eigenvalue
% (iii) calls qralg with smaller matrix to get another eigenvalue
% (iv) continue for all eigenvalues


A = hilb(4);
%A = diag(15:-1:1) + ones(15, 15); % part e

% call tridiag
T = tridiag(A);

[m,n] = size(T);
tvals = zeros(m-1,1);
eigens = zeros(m,1);
errvals = [];
for i=1:m-1
    
    % store T(m,m-1)
    tvals(i) = abs(T(end,end-1));
    
    % get one eigenvalue with qralg
    [Tnew, errvalsi] = qralg(T);
    
    % track the error value terms
    errvals = [errvals; errvalsi];
    
    % hold the eigenvalues
    eigens(i) = Tnew(end,end);
    
    % deflate matrix to get next eig
    T = Tnew(1:end-1, 1:end-1);
    
end
% grab the last eigenvalue
eigens(m) = T(1,1);
% flip to go from greatest to least
eigens = flipud(eigens);
eigens

% plot the error
semilogy(errvals)
xlabel('QR iteration')
ylabel('|T_{m,m-1}|')
title('Eigenvalue decomposition convergence')