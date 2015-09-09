% problem 5 

dx = 0.1;
dt = 0.01;
tend = 1;
xend = 1;

nnodes = 1/dx + xend;


% mass matrix
Mij = gallery('tridiag', nnodes, dx/6 - dt/(2*dx), 2*dx/3 + dt/dx, dx/6 - dt/(2*dx));

% B matrix
Bij = gallery('tridiag', nnodes, dx/6 + dt/(2*dx), 2*dx/3 - dt/dx, dx/6 + dt/(2*dx));

% create the initial condition
xnodes = dx*[0:nnodes-1]';
alpha0 = 0.5 - abs(xnodes-0.5);

% enforce the boundary conditions
Bij(1,1) = 1;
Bij(1,2) = 0;
Bij(end,end) = 1;
Bij(end, end-1) = 0;
Mij(1,1) = 1;
Mij(1,2) = 0;
Mij(end,end) = 1;
Mij(end, end-1) = 0;


figure()
plot(xnodes, alpha0);

tcur = 0;
alpha_old = alpha0;
while (tcur<tend)
    tcur = tcur + dt;
    
    intmat = Mij\Bij;
    alpha = intmat*alpha_old;
    disp(tcur)
    
    alpha_old = alpha;
    
end

figure()
plot(xnodes, alpha);
xlabel('x')
ylabel('u')
%ylim([0,0.5])

max(alpha)

    
    