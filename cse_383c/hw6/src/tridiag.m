function T = tridiag(A)
%
% D. Pederson
%
% Reduce a real symmetric m x m matrix to tridiagonal form by orthogonal
% similarity transformation

[m,n] = size(A);

if (m ~= n)
    error('Matrix must be square!');
end

for k=1:m-2
    x = A(k+1:m, k);
    v = zeros(size(x));
    v(1) = sign(x(1))*norm(x);
    v = v+x;
    v = v/norm(v);
    
    A(k+1:m, k:m) = A(k+1:m, k:m) - 2*v*(v'*A(k+1:m, k:m));
    A(1:m, k+1:m) = A(1:m, k+1:m) - 2*(A(1:m, k+1:m)*v)*v';
    
end

T = A;