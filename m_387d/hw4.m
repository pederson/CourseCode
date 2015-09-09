% problem 6

dt = 0.006 ;
dx = 0.1;
xstart = 0;
xend = 1;
tstart = 0;
tend = 1;

t = tstart:dt:tend;
x = xstart:dx:xend;

nnodes = 1/dx + xend;

% create the initial condition
xnodes = dx*[0:nnodes-1]';
u0 = 0.5 - abs(xnodes-0.5);

figure()
plot(xnodes, u0);

tcur = tstart;
uold = u0;
ucur = uold;
while (tcur<tend)
    tcur = tcur + dt;
    
    for i=2:nnodes-1
        ucur(i) = uold(i) + dt/(dx*dx)*(uold(i+1) - 2*uold(i) + uold(i-1));
    end
    
    disp(tcur)
    uold = ucur;
    
end

figure()
plot(xnodes, ucur);
xlabel('x')
ylabel('u')
%ylim([0,0.5])

max(ucur)
