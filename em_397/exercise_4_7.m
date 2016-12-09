% Exercise 4.7
% Bayesian Inverse Problems
% D. Pederson

a = -1;
D = sqrt(2*pi);
nsamples = 10000;
pim = @(x) heaviside((x-a)).*exp(-x.^2/2);
q = @(x) 1/sqrt(2*pi)*exp(-x.^2/2);

% Rejection-Acceptance sampling
msamp = zeros(nsamples,1);
n=0;
while n<nsamples
    
    % draw from proposal density
    m = randn(1);
    
    % compute acceptance probability
    alpha = pim(m)/(D*q(m));
    
    % draw from uniform dist
    coin = rand(1);
    
    if alpha > coin
        n=n+1; % accept
        msamp(n) = m;
    end
end


histogram(msamp)