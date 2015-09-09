%% testing_hw4-1
% Dylan Pederson
clear; clc;

for i=2:5
    
    N = 2^i;
    
    a = steady_diffusion_advection_1d(@c_ex, @f_ex, @v_ex, @p_ex, phi_ex(0), phi_ex(1), N);
    h = 1/(N-1);
    
    
    % since phi_approximate(grid points) = a, the coefficients, I can compute the error at those points.
    x = [0:h:1]';

    % evaluate analytic expression for phi(x) at the grid points. 
    phiex = phi_ex(x); 

    %Visualize solution
    err = norm(phiex - a)/norm(phiex);
    fprintf('[%d] points: relative error is %e\n', N, err);

    clf
    hold on;
    plot(x,phiex,'r', 'LineWidth',10); 
    plot(x,a,'b.','LineWidth',10);
    grid on; 
    title('Red line: exact solution, Blue dots: approximate','FontSize',20);

    hold off;
    
    pause

end