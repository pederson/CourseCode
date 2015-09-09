% hw2 problem 5

h = 0.2;
u0 = 1;
u1 = 1;

vstenc = [-1, -h/2, 1, -h/2];
ustenc = [-h/2, -1, -h/2, 1];
mcol = 2;

sz = 2*(1/h+1);
A = zeros(sz);


% apply stencil
for i=3:sz-1
    if mod(i,2)==1
        A(i,i-2:i+1) = vstenc;
    else
        A(i,i-3:i) = ustenc;
    end
end
A(1,1:2) = [ 1, -h/2];

size(A) 

% insert boundary conditions in the matrix
A(2,:) = zeros(1,sz);
A(2,2) = 1;
A(end,:) = zeros(1,sz);
A(end,end) = 1;

A

% insert boundary conditions in the rhs
rhs = zeros(sz,1);
rhs(2) = 1;
rhs(end) = 1;

x = A\rhs


% richardson extrapolation part
h_ = 0.2;
yhover2 = [1, 1.03, 1.102, 1.217, 1.382, 1.0]';
yh = [1, 1.04, 1.123, 1.25, 1.43, 1]';

% form rhs

rhs2 = zeros(length(yh)*2,1);
for i=1:length(yh)
    rhs2(2*i-1) = yh(i);
    rhs2(2*i) = yhover2(i);
end
rhs2

% form the matrix 
B = zeros(length(yh)*2);
stenc = [1, h_^2; 1, h_^2/4];
for i=1:length(yh)
    B(2*i-1:2*i, 2*i-1:2*i) = stenc;
end

size(B)
B
refined = B\rhs2;
refined
