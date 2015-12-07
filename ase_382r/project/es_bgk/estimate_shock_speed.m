function [] = estimate_shock_speed(txtfile)

[x,t0,t05,t10,t15,t20,t25,t30,t35,t40,t45,t50] = import_profile(txtfile);
obj = dlmread(txtfile,',', 3, 0);
x = obj(:,1);

alltime = [t0,t05,t10,t15,t20,t25,t30,t35,t40,t45,t50];
% find position where properties midway between max and min
pmax = max(t50);
pmin = min(t50);
p50 = (pmax+pmin)/2;

[rows, cols] = size(alltime);
pos50 = zeros(1,cols);
t = zeros(1,cols);

for i=2:cols
    pos50(i) = x(abs(alltime(:,i)-p50) == min(abs(alltime(:,i)-p50)));
    t(i) = (i-1)*5;
end

figure()
plot(pos50,t,'k-o')

sp = diff(pos50)./diff(t);
speed = mean(sp);
disp(['estimated shock speed is ',num2str(speed)])
disp(['variance of shock speed is ',num2str(var(sp))])

end