function a = steady_diffusion_advection_1d( c, f, v, p, phi0, phi1, N)
% Dylan Pederson
%
% 1D example for solving diffusion equation using finite elements (linear piecewise polynomials
% as a basis function. 
% solve - d/dx [ c(x) d/dx phi(x) ] + v(x) d/dx [phi(x)] + p(x) phi(x) = f(x), phi(x=0)=phi0, phi(x=1)=phi1   % boundary
% conditions. 

disp('hi')
h = 1/(N-1);  % interval (x_i+1 - x_i)  or element length
x = [0:h:1]'; % grid points - they define the intervals for the piecewise functions. 

% basis functions in [0,h]
b0 = @(s) 1 - s/h;
b1 = @(s) s/h;


% their derivatives
db0 = -1/h;  
db1 =  1/h;

% initialize memory for matrix and right hand side (rhs)
K= zeros(N);
rhs= zeros(N,1);


% Loop over over the intervals and compute integrals for rhs and K. 
for j=1:N-1 
  
  s0 = x(j);  % starting point for interval j. 
  
  rhs0 = quadl(@(s) f(s0+s).* b0(s), 0, h); 
	rhs1 = quadl(@(s) f(s0+s).* b1(s), 0, h);
	rhs(j)   = rhs(j)   + rhs0;
	rhs(j+1) = rhs(j+1) + rhs1;

  % diagonal terms
  K(j,j)     = K(j,j)     + quadl(@(s) db0*db0*c(s0+s) , 0, h);
  myfun = @(s) b0(s).*db0.*v(s0+s);
  K(j,j)     = K(j,j)     + quadl(myfun, 0, h);
  myfun = @(s) b0(s).*b0(s).*p(s0+s)';
  K(j,j)     = K(j,j)     + quadl(myfun , 0, h);
  
  K(j+1,j+1) = K(j+1,j+1) + quadl( @(s) db1*db1*c(s0+s) , 0,h );
  myfun = @(s) b1(s).*db1.*v(s0+s);
  K(j+1,j+1) = K(j+1,j+1) + quadl(myfun, 0, h);
  myfun = @(s) b1(s).*b1(s).*p(s0+s)';
  K(j+1,j+1) = K(j+1,j+1) + quadl(myfun, 0, h);
  
  % off diagonal terms
  tmp = quadl( @(s) db0*db1*c(s0+s) , 0,h );
  myfun = @(s) b0(s).*db1.*v(s0+s);
  tmp = tmp + quadl(myfun, 0, h);
  myfun = @(s) + b0(s).*b1(s).*p(s0+s)';
  tmp = tmp + quadl(myfun, 0, h);
  K(j,j+1)   = K(j,j+1)   + tmp;
  
  
  tmp = quadl( @(s) db1*db0*c(s0+s) , 0,h );
  myfun = @(s) b1(s).*db0.*v(s0+s);
  tmp = tmp + quadl(myfun, 0, h);
  myfun = @(s) + b1(s).*b0(s).*p(s0+s)';
  tmp = tmp + quadl(myfun, 0, h);
  K(j+1,j)   = K(j+1,j)   + tmp;
  
 
end

% Here we impose the boundary conditions. Since our solution phi(x) has to satisfy the boundary
rhs = rhs - K(:,1)*phi0 - K(:,end)*phi1;  % modify the right hand side appropriately
K = K(2:end-1,2:end-1);                   % solve only for the coefficients of the basis
                                          % functions i=2,...,N-1

rhs = rhs(2:end-1);   % select the correct rhs

a = K\rhs;            % solve the linear system
a = [phi0;a;phi1];    % Collect all the coefficients for the basis functions.

end
