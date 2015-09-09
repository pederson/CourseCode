%% testing_hw4-2
% Dylan Pederson
clear; clc;

for i=2:5
    
    N = 2^i;
    
    % find the a coefficients for N points
    a = steady_diffusion_advection_1d(@c_ex, @f_ex, @v_ex, @p_ex, phi_ex(0), phi_ex(1), N);
    
    % calculate phi on the 2N grid
    h = 1/(2*N-1);
    phi_2N = zeros(2*N, 1);
    for j=1:2*N
        ind = round(j/2);
        if mod(j,2) == 1 || j == 2*N
            phi_2N(j) =  a(ind);
        else
            phi_2N(j) = a(ind)/2 + a(ind+1)/2;
        end
    end
    
    % find the a coefficients for 2N points
    a2N = steady_diffusion_advection_1d(@c_ex, @f_ex, @v_ex, @p_ex, phi_ex(0), phi_ex(1), 2*N);
    
    % compare residual
    rhoN = max(abs(phi_2N - a2N));
    err = norm(phi_2N - a2N)/norm(phi_2N);
    fprintf('[%d] points: rhoN is %e\t relative error is %e\n', N, rhoN, err);

end