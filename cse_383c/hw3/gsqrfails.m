clear all; clear globals; 
m=500; n=100;

% Test different  type of matrices
Ag=randn(m,n);  % a random Gaussian matrix.
Ao = ones(m,m)+1e-11*rand(m,m); % a nearly rank-deficient matrix.
conditionnumber = 1e8;
As = gallery('randsvd',[m,n],conditionnumber); % a random ill-conditioned matrix.

% Select matrix by editing the line below: A = one of Ag, Ao, As
% For Ag it works perfectly fine. But for Ao or As it fails.
A=As;
A=A(:,1:n);

% THE GRAM-SCHMIDT ALGORITHM BEGINS HERE;
q = A(:,1); 
q = q/norm(q);
Q = [q];
for j=2:n
  aj = A(:,j);
  q = aj - Q*(Q'*aj);   % q=(I-Q*Q')aj
  q = q/norm(q);
  Q = [Q, q];
end
% AND ENDS HERE.

% Checking errors

fprintf('\n----------------------------------------\n');
fprintf('The condition number of the test matrix is %1.2g\n', cond(A));
fprintf('Test accuracy of our simple Gram-Schmidt orthogonalization implementation\n');

fprintf('Orthogonality error: %1.2e \t (should be O(1E-15)) \n', norm(Q'*Q-eye(n)));
z = randn(m,1);
z = z-Q*(Q'*z);
fprintf('  Range space error: %1.2e \t (should be O(1E-15)) \n', norm(A'*z));

[Q,~,~]=qr(A,0);
fprintf('\nErrors using MATLAB''s stable orthogonalization algorithm\n');
fprintf('Orthognoality error: %1.2e\n', norm(Q'*Q-eye(n)));
z = randn(m,1);
z = z-Q*(Q'*z);
fprintf('  Range space error: %1.2e\n', norm(A'*z));



%Exercise: what is the R factor and how you further check the error?