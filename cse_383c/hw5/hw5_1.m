function [] = hw5_1(A,b, epsilon, xstar)

[m,n] = size(A);
rank = min([m,n]);

%% part A
% compute the reduced SVD and plot singular values
[U, S, V] = svd(A, 0);
sing_vals = diag(S);
semilogy(sing_vals);
title('Singular values of A')
xlabel('n')
ylabel('singular value')

%% part B
% explain whether A is rank deficient (or nearly)


%% part C
% test svd, qr w/ pivot, qr w/o pivot, and rangeA.m function
% to find the best rank-r approximation to Range(A)

% QR w/ pivot
[Qpiv, Rpiv, Ppiv] = qr(A, 0);

% QR w/o pivot
[Q,R] = qr(A);

% rangeA
%QrangeA =rangeA(  @(x)A*x,  @(x)A'*x,  m, n, n);

for r = n*[1/16, 1/8, 1/4, 1/2]
    disp(['r: ',num2str(r)])
    
    % svd
    Ur = U(:,1:r);
    Sr = S(1:r, 1:r);
    Vr = V(:,1:r);

    % qr w/ pivot
    Qpivr = Qpiv(:,1:r);
    Rpivr = Rpiv(1:r, :);

    % qr w/o pivot
    Qr = Q(:,1:r);
    Rr = R(:,1:r);

    % rangeA.m
    QrangeA =rangeA(  @(x)A*x,  @(x)A'*x,  m, n, r);
    
    % error analysis
    disp(['SVD error: ',num2str(norm(A - Ur*Sr*Vr'))])   % SVD
    disp(['QR w/pivot error: ',num2str(norm(A - Qpivr*Rpivr))])   % QR w/pivot
    disp(['QR w/o pivot error: ',num2str(100)])   % QR w/o pivot
    disp(['rangeA error: ',num2str(100)])   % rangeA

end

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
    
    % relative cond # A
    
    
    % relative cond # b
    
    
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
plot(r, relcond_A, 'k')
title('relative condition number (pertub A)')
xlabel('Rank r')
ylabel('relative cond. number')

% (ii)
% plot relative cond, pertub b
figure()
plot(r, relcond_b, 'k')
title('relative condition number (pertub b)')
xlabel('Rank r')
ylabel('relative cond. number')

% (iii)
% plot the norm |xr|/|xstar| as fn(r)
figure()
semilogy(r, normx, 'k')
title('scaled norm of x_r')
xlabel('Rank r')
ylabel('|x_r|/|x_{*}|')

% (iv)
% plot residual as fn(r)
figure()
semilogy(r, resid, 'k')
title('residual')
xlabel('Rank r')
ylabel('residual')

% (v)
% plot the error as fn(r)
figure()
semilogy(r, error, 'k')
title('error in x_r')
xlabel('Rank r')
ylabel('error')

% (vi)
% plot the residual v error for each r using loglog
figure()
loglog(error, resid, 'k')
title('residual v. error')
xlabel('error')
ylabel('residual')

%% part E
% repeat part D but with regularized svd (beta = sigma_r)

