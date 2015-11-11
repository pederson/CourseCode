function Q=rangeA(A,At,m,n,r,q)
% function Q=rangeA(A,At,m,n,r,q)
%
%  A: function handle for A*x
% At: function handle for A'*x
%  m: size(A,1)
%  n: size(A,2)
%  r: target rank
%  q: optional value for power iteration (leave it unspecified)
%
%
% EXAMPLE
%  m=200;  n=100; r=50;
%  A =randn(m,n);
%  Q =rangeA(  @(x)A*x,  @(x)A'*x,  m, n, r);
%  norm(A-Q*(Q'*A))/norm(A);


if nargin<1, selfcheck; Q=[]; return; end;
if nargin<6, q=3; end;

G = randn(n,r+10);
Y = A(G);
[Q,~] = qr(Y,0);

% power iteration correction
for j=2:q
  Y=A(At(Q));
  [Q,~] = qr(Y,0);
end

% output
Q=Q(:,1:r);


function selfcheck
m=100; r=50; 
A = gallery('randsvd',m,1e10);
A=A(:,1:end-1); n=m-1;
[U,S,~]=svd(A); U=U(:,1:r);
Q =rangeA(@(x)A*x,@(x)A'*x,m,n,r);




