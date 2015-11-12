%% CSE 383C Homework 5 - Dylan Pederson

%% Question 1: Least-Squares Problems

% loop over each problem set
for p=1:5
    % skip part d and e for dataset 1 and 5
    if p==1 || p==5
        do_part_d = false;
        do_part_e = false;
    else
        do_part_d = true;
        do_part_e = true;
    end
    
    % label for plots
    problemstring = ['Problem ',num2str(p)];
    
    % load dataset
    [A, b, epsilon, xstar] = problem(p);
    
    % matrix properties
    [m,n] = size(A);
    rank = min([m,n]);

    %% part A
    % compute the reduced SVD and plot singular values
    [U, S, V] = svd(A, 0);
    sing_vals = diag(S);
    
    % plot singular values
    figure()
    semilogy(sing_vals, 'b');
    hold on
    semilogy([1, length(sing_vals)], [epsilon, epsilon], 'k--')
    title(['Singular values A: ', problemstring])
    xlabel('n')
    ylabel('singular value')
    legend('Singular values','epsilon')
    
    
    %% part B
    % explain whether A is rank deficient (or nearly).
    % If any singular values of A are smaller than epsilon
    % then the matrix is rank deficient
    if (nnz(sing_vals < epsilon) > 0)
        disp([problemstring,': Matrix is rank-deficient or nearly rank-deficient'])
    else
        disp([problemstring,': Matrix is NOT rank-deficient'])
    end

    %% part C
    % test SVD, QR w/ pivot, QR w/o pivot, and rangeA.m function
    % to find the best rank-r approximation to Range(A)
    % 
    % Here I am defining the error as:
    %
    % $$ error_{Range(A)} = ||A - Q_r Q_r' A||
    %
    % If we have a matrix A and a basis for Range(A) called Q
    % Then A z brings any vector z into Range(A).
    % A complementary projector I-QQ' brings a vector to Null(A)
    % If we consecutively apply A and the complementary projector I-QQ' to
    % a vector z then
    % (I-QQ')Az = 0
    % because if a vector in Range(A) is projected onto Null(A) then the
    % result must be the zero vector.
    % Taking norms, we find
    %
    % $$  ||(A - Q Q' A)z|| = 0
    %
    % the LHS is just the matrix norm so
    %
    % $$  ||A - Q Q' A|| = 0
    %
    % We can measure how "good" a basis is with this norm 
    
    % obtain the full factorizations first (where possible)
    % QR w/ pivot
    [Qpiv, ~] = qr(A, 0);

    % QR w/o pivot
    [Q,~] = qr(A);
    
    % initialize stuff to hold data
    svderr = zeros(4,1);
    qrpiverr = zeros(4,1);
    qrerr = zeros(4,1);
    rangeAerr = zeros(4,1);
    rs = n*[1/16, 1/8, 1/4, 1/2];
    for i=1:length(rs)
        r = rs(i);

        % truncate svd
        Ur = U(:,1:r);

        % truncate qr w/ pivot
        Qpivr = Qpiv(:,1:r);

        % truncate qr w/o pivot
        Qr = Q(:,1:r);

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
    xlabel('Rank r')
    ylabel('error in Range(A)')
    
    
    

    
    %% part D
    % use truncated SVD to solve least squares
    %
    % Discussion:
    %
    %     Problem 2: 
    %          - The relative condition numbers trend smoothly upwards as
    %          does the norm of x_r
    %          - The residual decreases and plateaus for a while around
    %          r=200-400 then decreases more
    %          - The error in x_r has a minimum point around r=200
    %          - Seems that the optimal choice is r~250
    %
    %     Problem 3:
    %          - Sharp discontinuities in relative condition number around
    %          r=100. This makes sense given the sharp drop in singular
    %          value around r=100
    %          - Sharp discontinuities in most other values after r=100
    %          leading to extreme errors
    %          - A good choice of r is r<100 (perhaps r=95)
    %
    %     Problem 4:
    %          - Same sharp discontinuity as in (3), but this time at r~75 
    %          - A good choice would be r<75
    
    if (do_part_d)
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

        % (ii)
        % plot relative cond, pertub b
        figure()
        semilogy(r, relcond_b, 'k+')
        title(['Rel. Cond. # (pertub b) (Truncated): ',problemstring])
        xlabel('Rank r')
        ylabel('relative cond. number b')

        % (iii)
        % plot the norm |xr|/|xstar| as fn(r)
        figure()
        semilogy(r, normx, 'k*')
        title(['Norm of x_r (Truncated): ',problemstring])
        xlabel('Rank r')
        ylabel('||x_r||/||x_{*}||')

        % (iv)
        % plot residual as fn(r)
        figure()
        semilogy(r, resid, 'ko')
        title(['Residual (Truncated): ',problemstring])
        xlabel('Rank r')
        ylabel('||b-Ax_r||/||b||')

        % (v)
        % plot the error as fn(r)
        figure()
        semilogy(r, error, 'kx')
        title(['Error in x_r (Truncated): ',problemstring])
        xlabel('Rank r')
        ylabel('||x_r-x_{star}||/||x_{star}||')

        % (vi)
        % plot the residual v error for each r using loglog
        figure()
        loglog(error, resid, 'kd')
        title(['residual v. error (Truncated): ',problemstring])
        xlabel('||x_r-x_{star}||/||x_{star}||')
        ylabel('||b-Ax_r||/||b||')
    end
    
    
    
    %% part E
    % repeat part D but with regularized svd (beta = sigma_r)
    %
    % The regularized SVD solution is 
    %
    % $$ x = V(\beta I + \Sigma ^2)^{-1}\Sigma U^* b
    %
    % And this is truncated in the same way as the normal SVD
    %
    % Discussion:
    %
    %     Problem 2: 
    %          - The relative condition numbers trend upwards
    %          - The norm of x_r increases then plateaus around r=100. This
    %          makes sense because of the regularization condition
    %          - The residual and error in x_r steadily decrease
    %          - A good choice of r would come from a tradeoff b/w desired
    %          accuracy and speed
    %
    %     Problem 3:
    %          - Sharp discontinuities in relative condition number around
    %          r=100.
    %          - Sharp discontinuities in most other values after r=100
    %          - Residual is fairly constant. Despite large condition
    %          number, error in x_r is minimized just above r=100
    %          - A good choice of r is r>100 (perhaps r=105)
    %
    %     Problem 4:
    %          - Same sharp discontinuity as in part (d)
    %          - As noted before, norm of x_r flattens out for higher r due
    %          to the regularization. This causes other values to flatten
    %          out for higher r
    
    if do_part_e
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
    end
    
    
end

%% Question 2: Matrix-Free
% Algorithm:
%
%       We have an orthonormal basis Z (mxr) given from rangeA.m.
%       Let B = Z'A (rxn).
%       Then B' = A'Z, and the colums of B' (b'_j) can be obtained by
%       multiplying A'*z_j.
%       B can be reconstructed by transposing B'.
%       Then perform a QR factorization on B (B = WT).
%       Note that B = Z'A = WT, so then A = ZWT where Z and W are
%       orthonormal matrices, and T is upper triangular.
%       Then we can solve for x.
%       Ax = b, then ZWTx=b, so x = T^{-1}W'Z'b 
% 
% Selecting the best r:
%
%       In order to get the best r for this problem, we look at the
%       residual |b-Ax|/|b|. Because we have no other information (besides
%       the relative condition nubmers), we must choose the rank where the
%       residual is no longer decreasing (as much). From the plots we see
%       that the conditions numbers increase monotonically and the residual
%       decreases quickly, and comes to a steady value at a rank r~20. So
%       we want to reach that steady low residual value with as small of a
%       condition number as possible.
%       From the plots we see that r~20 satisfies these criteria

% load the matrix
o = matrixfree;
m = o.m;
n = o.n;
fullrank = min([m,n]);

allranks = 1:1:100;%fullrank;
relcond_A = zeros(size(allranks));
relcond_b = zeros(size(allranks));
resid = zeros(size(allranks));
for i = 1:length(allranks)
    targetrank = allranks(i);
    
    % form orthonormal basis
    Z = rangeA(@(x)o.A(x), @(x)o.At(x), m, n, targetrank);

    % form B = Z'A
    B = zeros(targetrank, n);
    for j=1:targetrank
        zj = Z(:,j);
        bj = o.At(zj);
        bj = bj';
        B(j, :) = bj;
    end

    % form the QR factorization of B
    [W,T] = qr(B);

    % now the decomp. is Z T' W'
    % the decomp. is Z W T
    % solve least squares system
    y = Z'*o.b;
    q = W'*y;
    x = T\q;
    
    % do some error analysis
    % condition number parameters
    kappa = cond(T);
    theta = acos(norm(Z*Z'*o.b)/norm(o.b));
    eta = norm(T)*norm(x)/norm(o.A(x));
    

    % relative cond # A
    relcond_A(i) = kappa/(eta*cos(theta));

    % relative cond # b
    relcond_b(i) = kappa + kappa^2*tan(theta)/eta;

    % residual
    resid(i) = norm(o.b - o.A(x))/norm(o.b);
end

% (i)
% plot relative cond, pertub A
figure()
semilogy(allranks, relcond_A, 'k.')
title(['Rel. Cond. # (pertub A) (Matrix-Free): '])
xlabel('Rank r')
ylabel('relative cond. number B')

% (ii)
% plot relative cond, pertub b
figure()
semilogy(allranks, relcond_b, 'k+')
title(['Rel. Cond. # (pertub b) (Matrix-Free): '])
xlabel('Rank r')
ylabel('relative cond. number b')

% (iv)
% plot residual as fn(r)
figure()
semilogy(allranks, resid, 'ko')
title(['Residual (Matrix-Free): '])
xlabel('Rank r')
ylabel('||b-Ax_r||/||b||')
