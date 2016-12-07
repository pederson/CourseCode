% Problem 4 - Gradient Descent
%
% D. Pederson
% EE 381K

% the dataset
n = 100;
A = randn(n,n);
[U,R] = qr(A);

Sigma1 = diag(ones(n,1));
Sigma2 = diag([ones(n/2,1);0.01*ones(n/2,1)]);
Sigma3 = diag([100*ones(10,1);ones(n-10,1)]);
X1 = U*Sigma1*U';
X2 = U*Sigma2*U';
X3 = U*Sigma3*U';
b1 = ones(n,1);
b2 = b1;
b3 = b1;

% gradient descent with constant step size
X = X2;

% choose gamma as near the largest eigenvalue
d = eigs(X);
gamma = d(1)*0.0003;
gamma = 1;


beta_0 = ones(n,1);

niters = 20;
fbeta = zeros(niters, 1);
beta = beta_0;
for i=1:niters
    fbeta(i) = 0.5*beta'*X*beta;
    beta = beta - gamma*(X*beta);
end

figure()
plot(1:niters, fbeta)
xlabel('n')
ylabel('f(\beta)')
title(['Gradient Descent (\gamma=',num2str(gamma),')'])
