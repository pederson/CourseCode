function [Tnew, errvals] = qralg(T)
%
% D. Pederson
%
% runs unshifted QR algorithm on a real tridiagonal matrix T

[m,n] = size(T);

errvals = zeros(m,1);
errvals(1) = abs(T(m,m-1));
ct=1;
while (abs(T(m,m-1)) >= 1e-12)

    % qr factorize T
    [Q, R] = qr(T);
    
    % update T
    T = R*Q;    
    
    % update counter
    ct = ct+1;
    
    % track the error term
    errvals(ct) = abs(T(m,m-1));
    
end

errvals = errvals(1:ct);
ct
Tnew = T;