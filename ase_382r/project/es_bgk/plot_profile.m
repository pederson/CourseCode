function [] = plot_profile(txtfile)

[x,t0,t05,t10,t15,t20,t25,t30,t35,t40,t45,t50] = import_profile(txtfile);

figure()
hold on
plot(x, t0, 'ks-')
plot(x, t10, 'bo-')
plot(x, t20, 'r^-')
plot(x, t30, 'm>-')
plot(x, t40, 'gx-')
plot(x, t50, 'k+-')
legend('t=0','t=10','t=20','t=30','t=40','t=50')
xlabel('x position')

end