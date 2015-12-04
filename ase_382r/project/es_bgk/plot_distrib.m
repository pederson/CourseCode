function [] = plot_distrib(txtfile)

[vx,vy,phi] = import_distribfile(txtfile);

npts = length(phi);
snpts = sqrt(npts);

vx = reshape(vx,[snpts,snpts]);
vy = reshape(vy,[snpts,snpts]);
phi = reshape(phi,[snpts,snpts]);

figure()
surf(vx,vy,phi)
xlabel('V_x/\eta_x')
ylabel('V_y/\eta_y')
zlabel('$$\hat{\phi}$$','Interpreter','Latex')

end