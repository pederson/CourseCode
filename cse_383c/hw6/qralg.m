function Tnew = qralg(T)
%
% D. Pederson
%
% runs unshifted QR algorithm on a real tridiagonal matrix T

[m,n] = size(T);

while (abs(T(m,m-1)) >= 1e-12)

    % qr factorize T
    [Q, R] = qr(T);
    
    % update T
    T = R*Q;    
    
end

Tnew = T;