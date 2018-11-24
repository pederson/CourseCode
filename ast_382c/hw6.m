%% Homework 6, Problem 2 (d)

kx = 0:0.05:(2*pi);

gamma = 5/3;
eps = 0.95;


ct = 0:0.01:1;

cs = (1+eps*cos(kx)).^((gamma-1)/2);

figure()
plot(kx, cs)
ylabel('c_s')

% v0 = -eps*(1+eps*cos(kx)).^((gamma-3)/2).*cos(kx);
v0 = 2/(gamma-1)*(1-cs);


% slope_cplus = (1+eps*cos(kx)).^((gamma-1)/2).*(1-eps*cos(kx).*(1+cos(kx)).^-1);
slope_cplus = v0 + cs;

slope_cminus = -slope_cplus;

figure()
plot(kx/pi, v0)

figure()
plot(kx/pi, 1./slope_cplus)


figure()
hold on
for i = 1:length(kx)
    plot([kx(i), kx(i)+pi/slope_cplus(i)]/pi, [0, pi], 'k-');
end
xlim([0, 2*pi]/pi)
set(gca, 'FontSize', 16)
xlabel('$ kx / \pi $', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$ k c_{s,0} t$', 'Interpreter', 'latex', 'FontSize', 20)