N = 1000; 
alpha = 0.1;
x = (1:N)/N;
t=1:100;

f0 = zeros(length(x));
for i=1:length(x)/2
    f0(i) = x(i)*2;
end
for i=length(x)/2+1:length(x)
    f0(i) = 2 - 2*x(i);
end


y = wave_1D_bounded(x, t, alpha, f0, N/10);

% figure()
% hold on
% for n=1:length(t)
%     plot(x, y(n,:))
% end

surf(x, t, y, 'EdgeColor','none')