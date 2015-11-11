function [] = hw5_1(A,b, epsilon, xstar, problemstring)
if nargin<5, problemstring = ''; end;
    
[m,n] = size(A);
rank = min([m,n]);

%% part A
% compute the reduced SVD and plot singular values
[U, S, V] = svd(A, 0);
sing_vals = diag(S);
figure()
hold on
semilogy(sing_vals, 'b.');
%semilogy([1, length(sing_vals)], [epsilon, epsilon], 'k--')
title(['Singular values A: ', problemstring])
xlabel('n')
ylabel('singular value')
%legend('Singular values','epsilon')

%% part B
% explain whether A is rank deficient (or nearly)
% if any singular values of A are smaller than epsilon
% then the matrix is rank deficient
if (nnz(sing_vals < epsilon) > 0)
    disp('Matrix is rank-deficient')
else
    disp('Matrix is NOT rank-deficient')
end

%% part C
% test svd, qr w/ pivot, qr w/o pivot, and rangeA.m function
% to find the best rank-r approximation to Range(A)

% QR w/ pivot
[Qpiv, Rpiv, Ppiv] = qr(A, 0);

% QR w/o pivot
[Q,R] = qr(A);

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
%     disp(['    SVD error: ',num2str(norm(A - Ur*Sr*Vr'))])   % SVD
%     disp(['    QR w/pivot error: ',num2str(norm(A - Qpivr*Rpivr))])   % QR w/pivot
%     disp(['    QR w/o pivot error: ',num2str(norm(I-Qr'*Qr))])   % QR w/o pivot
%     disp(['    rangeA error: ',num2str(100)])   % rangeA

end

r = n*[1/16, 1/8, 1/4, 1/2]';
figure()
hold on
plot(r, svderr, 'ko-')
plot(r, qrpiverr, 'b+-')
plot(r, qrerr, 'gd-')
plot(r, rangeAerr, 'r.-')
legend('SVD','QR w/pivot','QR w/o pivot','rangeA')

%% part D
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
plot(r, relcond_A, 'k.')
title('relative condition number (pertub A) (Truncated SVD)')
xlabel('Rank r')
ylabel('relative cond. number')

% (ii)
% plot relative cond, pertub b
figure()
plot(r, relcond_b, 'k+')
title('relative condition number (pertub b) (Truncated SVD)')
xlabel('Rank r')
ylabel('relative cond. number')

% (iii)
% plot the norm |xr|/|xstar| as fn(r)
figure()
semilogy(r, normx, 'k*')
title('scaled norm of x_r (Truncated SVD)')
xlabel('Rank r')
ylabel('|x_r|/|x_{*}|')

% (iv)
% plot residual as fn(r)
figure()
semilogy(r, resid, 'ko')
title('residual (Truncated SVD)')
xlabel('Rank r')
ylabel('residual')

% (v)
% plot the error as fn(r)
figure()
semilogy(r, error, 'kx')
title('error in x_r (Truncated SVD)')
xlabel('Rank r')
ylabel('error')

% (vi)
% plot the residual v error for each r using loglog
figure()
loglog(error, resid, 'kd')
title('residual v. error (Truncated SVD)')
xlabel('error')
ylabel('residual')

%% part E
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
plot(r, relcond_A, 'k.')
title('relative condition number (pertub A) (Regularized SVD)')
xlabel('Rank r')
ylabel('relative cond. number')

% (ii)
% plot relative cond, pertub b
figure()
plot(r, relcond_b, 'k+')
title('relative condition number (pertub b) (Regularized SVD)')
xlabel('Rank r')
ylabel('relative cond. number')

% (iii)
% plot the norm |xr|/|xstar| as fn(r)
figure()
semilogy(r, normx, 'k*')
title('scaled norm of x_r (Regularized SVD)')
xlabel('Rank r')
ylabel('|x_r|/|x_{*}|')

% (iv)
% plot residual as fn(r)
figure()
semilogy(r, resid, 'ko')
title('residual (Regularized SVD)')
xlabel('Rank r')
ylabel('residual')

% (v)
% plot the error as fn(r)
figure()
semilogy(r, error, 'kx')
title('error in x_r (Regularized SVD)')
xlabel('Rank r')
ylabel('error')

% (vi)
% plot the residual v error for each r using loglog
figure()
loglog(error, resid, 'kd')
title('residual v. error (Regularized SVD)')
xlabel('error')
ylabel('residual')
