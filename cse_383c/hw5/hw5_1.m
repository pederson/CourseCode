function [] = hw5_1(A,b, epsilon, xstar, problemstring)
if nargin<5, problemstring = ''; end;
    
[m,n] = size(A);
rank = min([m,n]);

%% part A
disp('Part A')
% compute the reduced SVD and plot singular values
[U, S, V] = svd(A, 0);
sing_vals = diag(S);
% figure()
semilogy(sing_vals, 'b');
hold on
semilogy([1, length(sing_vals)], [epsilon, epsilon], 'k--')
title(['Singular values A: ', problemstring])
xlabel('n')
ylabel('singular value')
legend('Singular values','epsilon')
%savefig(['p1_a_singvals_',problemstring,'.fig'])

%% part B
disp('Part B')
% explain whether A is rank deficient (or nearly)
% if any singular values of A are smaller than epsilon
% then the matrix is rank deficient
if (nnz(sing_vals < epsilon) > 0)
    disp('Matrix is rank-deficient')
else
    disp('Matrix is NOT rank-deficient')
end

%% part C
disp('Part C')
% test svd, qr w/ pivot, qr w/o pivot, and rangeA.m function
% to find the best rank-r approximation to Range(A)

% QR w/ pivot
[Qpiv, ~] = qr(A, 0);

% QR w/o pivot
[Q,~] = qr(A);

% rangeA
%QrangeA =rangeA(  @(x)A*x,  @(x)A'*x,  m, n, n);
svderr = zeros(4,1);
qrpiverr = zeros(4,1);
qrerr = zeros(4,1);
rangeAerr = zeros(4,1);
rs = n*[1/16, 1/8, 1/4, 1/2];
for i=1:length(rs)
    r = rs(i);
    disp(['r: ',num2str(r)])
    
    % svd
    Ur = U(:,1:r);
    %Sr = S(1:r, 1:r);
    %Vr = V(:,1:r);

    % qr w/ pivot
    Qpivr = Qpiv(:,1:r);
    %Rpivr = Rpiv(1:r, :);

    % qr w/o pivot
    Qr = Q(:,1:r);
    %Rr = R(:,1:r);

    % rangeA.m
    QrangeA =rangeA(  @(x)A*x,  @(x)A'*x,  m, n, r);
    
    % error analysis
    svderr(i) = norm(A - Ur*Ur'*A);
    qrpiverr(i) = norm(A - Qpivr*Qpivr'*A);
    qrerr(i) = norm(A - Qr*Qr'*A);
    rangeAerr(i) = norm(A - QrangeA*QrangeA'*A);
    
end

r = n*[1/16, 1/8, 1/4, 1/2]';
figure()
semilogy(r, svderr, 'ko-')
hold on
semilogy(r, qrpiverr, 'b+-')
semilogy(r, qrerr, 'gd-')
semilogy(r, rangeAerr, 'r.-')
legend('SVD','QR w/pivot','QR w/o pivot','rangeA')
title(['Low-rank approximations: ',problemstring])
%savefig(['p1_c_methods_',problemstring,'.fig'])

%% part D
disp('Part D')
% use truncated SVD to solve least squares
relcond_A = zeros(rank,1);
relcond_b = zeros(rank,1);
normx = zeros(rank,1);
resid = zeros(rank,1);
error = zeros(rank,1);
for r=1:rank
    Ur = U(:,1:r);
    Sr = S(1:r, 1:r);
    Vr = V(:,1:r);
    
    % solve the least squares problem
    xr = Vr*((Ur'*b)./diag(Sr));
    
    % condition number parameters
    kappa = Sr(1,1)/Sr(end,end);
    theta = acos(norm(Ur*Ur'*b)/norm(b));
    eta = Sr(1,1)*norm(xr)/norm(A*xr);
    
    % relative cond # A
    relcond_A(r) = kappa/(eta*cos(theta));
    
    % relative cond # b
    relcond_b(r) = kappa + kappa^2*tan(theta)/eta;
    
    % norm x
    normx(r) = norm(xr)/norm(xstar);
    
    % residual
    resid(r) = norm(b-A*xr)/norm(b);
    
    % error
    error(r) = norm(xr - xstar)/norm(xstar);
    
end

r=1:rank;
% (i)
% plot relative cond, pertub A
figure()
semilogy(r, relcond_A, 'k.')
title(['Rel. Cond. # (pertub A) (Truncated): ',problemstring])
xlabel('Rank r')
ylabel('relative cond. number B')
%savefig(['p1_d_relcond_A_',problemstring,'.fig'])

% (ii)
% plot relative cond, pertub b
figure()
semilogy(r, relcond_b, 'k+')
title(['Rel. Cond. # (pertub b) (Truncated): ',problemstring])
xlabel('Rank r')
ylabel('relative cond. number b')
%savefig(['p1_d_relcond_b_singvals_',problemstring,'.fig'])

% (iii)
% plot the norm |xr|/|xstar| as fn(r)
figure()
semilogy(r, normx, 'k*')
title(['Norm of x_r (Truncated): ',problemstring])
xlabel('Rank r')
ylabel('||x_r||/||x_{*}||')
%savefig(['p1_d_xr_',problemstring,'.fig'])

% (iv)
% plot residual as fn(r)
figure()
semilogy(r, resid, 'ko')
title(['Residual (Truncated): ',problemstring])
xlabel('Rank r')
ylabel('||b-Ax_r||/||b||')
%savefig(['p1_d_resid_',problemstring,'.fig'])

% (v)
% plot the error as fn(r)
figure()
semilogy(r, error, 'kx')
title(['Error in x_r (Truncated): ',problemstring])
xlabel('Rank r')
ylabel('||x_r-x_{star}||/||x_{star}||')
%savefig(['p1_a_singvals_',problemstring,'.fig'])

% (vi)
% plot the residual v error for each r using loglog
figure()
loglog(error, resid, 'kd')
title(['residual v. error (Truncated): ',problemstring])
xlabel('||x_r-x_{star}||/||x_{star}||')
ylabel('||b-Ax_r||/||b||')

%% part E
disp('Part E')
% repeat part D but with regularized svd (beta = sigma_r)
relcond_A = zeros(rank,1);
relcond_b = zeros(rank,1);
normx = zeros(rank,1);
resid = zeros(rank,1);
error = zeros(rank,1);
for r=1:rank
    Ur = U(:,1:r);
    Sr = S(1:r, 1:r);
    Vr = V(:,1:r);
    
    % regularization parameter
    beta = Sr(end,end);
    
    % solve the least squares problem
    extra = beta*eye(r) + Sr*Sr;
    extrainv = diag(1./diag(extra));
    xr = Vr*extrainv*Sr*Ur'*b;
    
    % condition number parameters
    kappa = Sr(1,1)/Sr(end,end);
    theta = acos(norm(Ur*Ur'*b)/norm(b));
    eta = Sr(1,1)*norm(xr)/norm(A*xr);
    
    % relative cond # A
    relcond_A(r) = kappa/(eta*cos(theta));
    
    % relative cond # b
    relcond_b(r) = kappa + kappa^2*tan(theta)/eta;
    
    % norm x
    normx(r) = norm(xr)/norm(xstar);
    
    % residual
    resid(r) = norm(b-A*xr)/norm(b);
    
    % error
    error(r) = norm(xr - xstar)/norm(xstar);
    
end

r=1:rank;
% (i)
% plot relative cond, pertub A
figure()
semilogy(r, relcond_A, 'k.')
title(['Rel. Cond. # (pertub A) (Regularized): ',problemstring])
xlabel('Rank r')
ylabel('relative cond. number B')

% (ii)
% plot relative cond, pertub b
figure()
semilogy(r, relcond_b, 'k+')
title(['Rel. Cond. # (pertub b) (Regularized): ',problemstring])
xlabel('Rank r')
ylabel('relative cond. number b')

% (iii)
% plot the norm |xr|/|xstar| as fn(r)
figure()
semilogy(r, normx, 'k*')
title(['Norm of x_r (Regularized): ',problemstring])
xlabel('Rank r')
ylabel('||x_r||/||x_{*}||')

% (iv)
% plot residual as fn(r)
figure()
semilogy(r, resid, 'ko')
title(['Residual (Regularized): ',problemstring])
xlabel('Rank r')
ylabel('||b-Ax_r||/||b||')

% (v)
% plot the error as fn(r)
figure()
semilogy(r, error, 'kx')
title(['Error in x_r (Regularized): ',problemstring])
xlabel('Rank r')
ylabel('||x_r-x_{star}||/||x_{star}||')

% (vi)
% plot the residual v error for each r using loglog
figure()
loglog(error, resid, 'kd')
title(['residual v. error (Regularized): ',problemstring])
xlabel('||x_r-x_{star}||/||x_{star}||')
ylabel('||b-Ax_r||/||b||')

