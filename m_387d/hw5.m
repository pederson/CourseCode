% hw 5 problem 6

dx = 0.02;
dt = 0.015;
tstart = 0;
tstop = 0.24;

t = tstart:dt:tstop;

% construct the initial solution
x = 0:dx:1;
for i=1:length(x)
    if x(i) > 0.25 && x(i) < 0.75
        u0(i) = 1;
    else
        u0(i) = -0.5;
    end
end

u = u0;
uold = u0;

u_alltime = zeros(length(x), length(t));

%%%%%% FDM method %%%%%%
%figure()
u_alltime(:,1) = u0;
for it = 2:length(t)
    for ix = 2:length(x)-1
        
        u(ix) = 0.5*(uold(ix+1) + uold(ix-1)) - dt/(2*dx)*(uold(ix+1)^2/2 - uold(ix-1)^2/2);
        u(1) = uold(end);
        u(end) = uold(1);
        
    end
    u_alltime(:,it) = u;
    uold = u;
    %plot(x, u, 'r');
end
u_FDM = u_alltime;


u = u0;
uold = u0;
u_alltime = zeros(length(x), length(t));

%%%%%% FVM method %%%%%%
%figure()
u_alltime(:,1) = u0;
for it = 2:length(t)
    for ix = 2:length(x)-1
        
        u(ix) = 0.5*(uold(ix+1) + uold(ix-1)) - dt*uold(ix)/(2*dx)*(uold(ix+1) - uold(ix-1));
        u(1) = uold(end);
        u(end) = uold(1);
        
    end
    u_alltime(:,it) = u;
    uold = u;
    %plot(x, u, 'b');
end
u_FVM = u_alltime;





u = u0;
uold = u0;
u_alltime = zeros(length(x), length(t));

%%%%%% DG method %%%%%%
%figure()
u_alltime(:,1) = u0;
for it = 2:length(t)
    for ix1 = 1:length(x)
        if u(ix1) > 0 
            fplus(ix1) = u(ix1)^2/2;
            fminus(ix1) = 0;
        else 
            fplus(ix1) = 0;
            fminus(ix1) = u(ix1)^2/2;
        end
    end
    for ix = 2:length(x)-1
        
        
        
        u(ix) =uold(ix) - dt/dx*(fminus(ix+1)-fminus(ix) + fplus(ix) - fplus(ix-1));
        u(1) = uold(end);
        u(end) = uold(1);
        
    end
    u_alltime(:,it) = u;
    uold = u;
    %plot(x, u, 'b');
end
u_DG = u_alltime;



figure()

for it = 1:length(t)
    plot(x, u_FDM(:,it), 'b');
    hold on;
    plot(x, u_FVM(:,it), 'r');
    plot(x, u_DG(:,it), 'g');
    hold off;
    legend('FDM', 'FVM', 'DG')
    title (['t = ',num2str(t(it))])
    pause;
end
    


        