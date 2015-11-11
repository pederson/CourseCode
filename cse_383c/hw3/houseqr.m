function [Q,R,U] = houseqr(A)

% Q: mxn matrix
% R: nxn matrix
% U: mxn matrix householder reflection vectors

    [m, n] = size(A);

    U = 0*A;
    R = A;
    
    for k=1:n
        w=R(k:m, k);
        %e1 = zeros(length(w), 1);
        %e1(1) = 1;
        %u = sign(w(1))*norm(w)*e1 + w;
        sigm = sign(w(1))*norm(w);
        u = w;
        u(1) = u(1) + sigm;
        u = u./norm(u);
        %b = 2/norm(u)^2;
        %b = 1/(sigm*u(1));
        
        R(k:m, k:n) = R(k:m, k:n) - 2*u*(u'*R(k:m, k:n));
        U(k:m, k) = u;
        
        R
    end

    % construct Q
    %Q = zeros(size(R));
    
    Q = A/R;
end
