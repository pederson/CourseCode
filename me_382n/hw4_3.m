%% testing_hw4-3
% Dylan Pederson
clear; clc;

bks = @(k,s) cos(k*s);
qx = @(x) 2*x -1; % change variables to go from [-1, 1]
qs = @(s) -cos(s); % transform to [0, pi]
xq = @(q) (q+1)/2;
sq = @(q) acos(-q);

for i=2:5
    
    N = 2^i;
    
    a = steady_diffusion_advection_cheb_1d(@c_ex, @f_ex, @v_ex, @p_ex, phi_ex(0), phi_ex(1), N);
    
    h = 1/(N-1);
    k=0:N-1;
    qj = qs(k*pi/(N-1));
    
    
    % since phi_approximate(grid points) = a, the coefficients, I can
    % compute the error at those points.
    phiapprox = zeros(N, 1);
    for j=1:N
            phiapprox(j) = sum(a'.*bks(k,sq(qj(j))));
    end

    % evaluate analytic expression for phi(x) at the grid points. 
    phiex = phi_ex(xq(qj))'; 

    %Visualize solution
    err = norm(phiex - phiapprox)/norm(phiex);
    fprintf('[%d] points: relative error is %e\n', N, err);

    clf
    hold on;
    plot(xq(qj),phiex,'r', 'LineWidth',10); 
    plot(xq(qj),phiapprox,'b.','LineWidth',10);
    grid on; 
    title('Red line: exact solution, Blue dots: approximate','FontSize',20);

    hold off;
    
    pause

end