function a = steady_diffusion_2d( c, f, g, N)
% Dylan Pederson
%
% 2D example for solving diffusion equation using finite elements (linear piecewise polynomials
% as a basis function on a regular grid.
% solve - div [ c(x,y) Grad phi(x,y) ] = f(x,y)
% x,y \in Omega
% Omega = [0,1]^2
% Gamma = boundary of Omega
%
% Boundary  conditions: phi(x,y) = g(x,y), point (x,y) on Gamma.
% Thermal conductivity, c(x,y)
% Heat source, f(x,y)
% The code will work with any f(x,y) and c(x,y). 
 
h = 1/(N-1);  
[x,y] = meshgrid(0:h:1); % grid points - they define the intervals for the piecewise functions. 

%% Define basis function on the reference domain

% 1D basis functions in [0,h]
b1d = { @(s) 1-s/h, @(s) s/h};

% combine 1D functions to build the 2D functions in the reference domain. 
% See Figure 3 in the lec13-fem.pdf for the exact expressions of the 2D
% functions. 
b = { 
  {@(x,y) b1d{1}(x).*b1d{1}(y),
   @(x,y) b1d{2}(x).*b1d{1}(y)},
  {@(x,y) b1d{1}(x).*b1d{2}(y),
   @(x,y) b1d{2}(x).*b1d{2}(y)}
};

%%
% their derivatives
db1d = [-1/h;1/h];

dbdx = {...
  {@(x,y) db1d(1)*b1d{1}(y),...
   @(x,y) db1d(2)*b1d{1}(y)},...
  {@(x,y) db1d(1)*b1d{2}(y),...
   @(x,y) db1d(2)*b1d{2}(y)}...
};
dbdy = {...
  {@(x,y) b1d{1}(x)*db1d(1),...
   @(x,y) b1d{2}(x)*db1d(1)},...
  {@(x,y) b1d{1}(x)*db1d(2),...
   @(x,y) b1d{2}(x)*db1d(2)}...
};

% initialize memory for matrix and right hand side (rhs)
K= zeros(N^2);
F= zeros(N^2,1);  
%%

% Convert i,j indices to  global indexing from 1,...,N^2
ij2g =@(i,j) (i-1)*N + j;

%g2j = @(k) mod(k-1,N)+1; g2ij =@(k)  [(k-g2j(k))/N+1, g2j(k)];

st = [0,1]; % auxiliary variable for looping over the vertices of the reference region.
%%
% Loop over over the regions (the finite elements) and compute integrals for rhs and K. 
for i=1:N-1,
  for j=1:N-1
  
    x0=x(i,j);  % x-coordinate of lower-left corner of region (i,j)
    y0=y(i,j);  % y-coordinate of lower-left corner of region (i,j)
    
    % loop over the four basis functions that are not zero in the i,j: they
    % are:  i,j; i+1,j; i,j+1; i+1,j+1;
    for mi=1:2     
      for mj=1:2
        
        % for each combination compute their global index
        m =  ij2g( i+st(mi), j+st(mj) );
        
        
        % COMPUTE  F (rhs vector)
        F( m) = F(m) + ...
          quad2d( @(x,y) f(x0+x,y0+y).*b{mi}{mj}(x,y),0,h,0,h);

        
        % COMPUTE MATRIX. 
        % for each m, we need to compute int( c grad b
        for ki=1:2
          for kj=1:2
            
            k = ij2g(i+st(ki),j+st(kj) );
            
            K(m,k) = K(m,k) + ...
              quad2d( ...
              @(x,y) c(x0+x,y0+y).*...
              (dbdx{mi}{mj}(x,y).*dbdx{ki}{kj}(x,y)+...
               dbdy{mi}{mj}(x,y).*dbdy{ki}{kj}(x,y)), ...
              0,h,0,h...
              );
          end
        end
        
        
      end
    end

  end
end
 
% boundary conditions
% first identify the boundary nodes. For that we will use the coordinates
% x,,y
boundary = (x==0 | x==1) | (y==0 | y==1);   % display boundary  to see what it looks like
interior =~boundary(:);  % make it a vector

% modify the rhs by subtracting out the a vectors
for i=1:N
    for j=1:N
        ind = ij2g(i,j);
    if boundary(ind)
        ind = ij2g(i,j);
        F = F - K(:,ind)*g(x(i,j),y(i,j));
    end
    end
end
K_interior=K(interior,interior); % pick all the interior nodes from K and F
F_interior = F(interior);
a_interior = K_interior\F_interior; % solve.
a = zeros*x; 
a=a(:);
a(interior)=a_interior;
a = reshape(a,N,N);
for i=1:N
    for j=1:N
        if i==1 || j==1 || i==N || j==N
            a(i,j) = g(x(i),y(i));
        end
    end
end

