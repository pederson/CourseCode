% ORTHOGONAL MATCHING PURSUIT
%
% D. Pederson
% EE 381K

k = 5;

X = X3;
y = y3;

Xtest = X3test;
ytest = y3test;

tic
I = [];
[rows,cols] = size(X);
yres = y;
for i=1:k
   
    % construct the new I
    maxcol = 0;
    maxval = 0;
    for j=1:cols
        val = yres'*X(:,j);
        if val > maxval
            maxval = val;
            maxcol = j;
        end
    end
    I = [I,maxcol];
    
    % construct the new Xi
    Xi = X(:,I);
    
    % solve the least squares problem
    bi_hat = Xi\y;
    
    % calculate the new residual
    yres = y - Xi*bi_hat;
    
    % repeat
end

b = zeros(cols, 1);
b(I) = bi_hat;
toc

% report the findings

% sparsity pattern
I

% regression error
norm(X*b-y,2)

% testing error
norm(Xtest*b-ytest,2)


