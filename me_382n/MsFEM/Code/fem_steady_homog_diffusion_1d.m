function [a] = fem_steady_homog_diffusion_1d(K, bc_left, bc_right)
% Dylan Pederson
% Miles Gray
%
% 1D diffusion equation using finite elements (linear piecewise polynomials
% as a basis function. 
% solve - d/dx [ c(x) d/dx phi(x) ] = 0, phi(x=0)=phi0, phi(x=1)=phi1   % boundary
% conditions. 
% The K matrix is known and passed in as a parameter

% Here we impose the boundary conditions. Since our solution phi(x) has to satisfy the boundary
rhs = -K(:,1)*bc_left -K(:,end)*bc_right;  % create the right hand side appropriately
K = K(2:end-1,2:end-1);                   % solve only for the coefficients of the basis
                                          % functions i=2,...,N-1

rhs = rhs(2:end-1);   % select the correct rhs

a = K\rhs;            % solve the linear system
a = [bc_left;a;bc_right];    % Collect all the coefficients for the basis functions.

end
