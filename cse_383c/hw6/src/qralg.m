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

    %{
    % uncomment this part for unshifted qr method
    % qr factorize T
    [Q, R] = qr(T);
    
    % update T
    T = R*Q;
    %}
    
    %
    % wilkinson shift method
    delta = (T(end-1, end-1) - T(end, end))/2;
    mu = T(end,end) - sign(delta)*T(end, end-1)^2/(abs(delta) + sqrt(delta^2 + T(end, end-1)^2));
    % qr factorize
    [Q, R] = qr(T - mu*eye(n));
    
    % recompose T
    T = R*Q + mu*eye(n);
    %}
    
    % update counter
    ct = ct+1;
    
    % track the error term
    errvals(ct) = abs(T(m,m-1));
    
end

errvals = errvals(1:ct);
Tnew = T;