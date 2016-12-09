% Exercise 4.8
% Bayesian Inverse Problems
% D. Pederson

sigma = 0.1;
delta = 1;
D = sqrt(2*pi);
nsamples = 10000;
pim = @(x) exp(-1/(2*sigma^2)*(sqrt(x(1).^2+x(2).^2)-1)^2 - 1/(2*delta^2).*(x(2).^2-1)^2);
q = @(x) 1/sqrt(2*pi)*exp(-x'*x/2);

% Rejection-Acceptance sampling
msamp = zeros(nsamples,2);
n=0;
while n<nsamples
    
    % draw from proposal density
    m = randn(2,1);
    
    % compute acceptance probability
    alpha = pim(m)/(D*q(m));
    
    % draw from uniform dist
    coin = rand(1);
    
    if alpha > coin
        n=n+1; % accept
        msamp(n,:) = m;
    end
end

x = -2:0.1:2;
y = x;
z = zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        z(i,j) = pim([x(i),y(j)]);
    end
end
contour(x,y,z)
hold on
for i=1:10:nsamples
    plot(msamp(i,1),msamp(i,2),'ro')
end