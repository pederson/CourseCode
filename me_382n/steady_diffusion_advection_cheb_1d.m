function a = steady_diffusion_advection_cheb_1d( c, f, v, p, phi0, phi1, N)
% Dylan Pederson
%
% 1D example for solving diffusion equation using chebyshev collocation 
% solve - d/dx [ c(x) d/dx phi(x) ] + v(x) d/dx [phi(x)] + p(x) phi(x) = f(x), phi(x=0)=phi0, phi(x=1)=phi1   % boundary
% conditions. defined on [0,1]

h = 1/(N-1);
k=0:N-1;

qx = @(x) 2*x -1; % change variables to go from [-1, 1]
qs = @(s) -cos(s); % transform to [0, pi]
xq = @(q) (q+1)/2; 

qj = qs(k*pi/(N-1));
sq = @(q) acos(-q);

% basis functions yay
bk = @(k, s) cos(k*s);
dbk = @(k, s) -k.*sin(k*s);
dbk2 = @(k, s) -k.*k.*cos(k*s);

% derivative of c
c_prime = cheb1_dfds(c(0:h:1))';
% adjust for mapping weight
c_prime = 2*c_prime;

%{
figure()
hold on
plot(0:h/5:1,c(0:h/5:1),'r')
plot(0:h:1,c_prime,'b')

pause
%}

% initialize memory for matrix
K = zeros(N);

% construct the K matrix elements
for i=1:N
    for j=1:N
        if i==1 || i==N % for Boundary Conditions
            K(i,j) = bk(k(j), sq(qj(i)));
            continue;
        end
        
        %-c(xq(qj(i)))*4*-qj(i)/(1-qj(i)*qj(i))^(3/2)*dbk(k(j), sq(qj(i))) ...
        K(i,j) = -c_prime(i)*2/sqrt(1-qj(i)*qj(i))*dbk(k(j), sq(qj(i))) ...
                 -c(xq(qj(i)))*4/(1-qj(i)*qj(i))*dbk2(k(j), sq(qj(i))) ...
                 -c(xq(qj(i)))*4*(-qj(i))/(1-qj(i)*qj(i))^(3/2)*dbk(k(j), sq(qj(i))) ...
                 +v(xq(qj(i)))*2/sqrt(1-qj(i)*qj(i))*dbk(k(j), sq(qj(i))) ...
                 +p(xq(qj(i)))*bk(k(j), sq(qj(i)));
    end
end

% construct the rhs vector
rhs = f(xq(qj))';

% Here we impose the boundary conditions. Since our solution phi(x) has to
% satisfy the boundary
rhs(1) = phi0;
rhs(N) = phi1;

a = K\rhs;            % solve the linear system

end
