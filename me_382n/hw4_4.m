%% testing_hw4-4
% Dylan Pederson
clear; clc;

for i=2:5
    
    N = 2^i;
    
    a = steady_diffusion_2d(@c_ex_2d, @f_ex_2d, @g_ex_2d, N);
    h = 1/(N-1);

    % vizyoolize
    clf; surf(a),shading interp,colormap default;

    pause

end

