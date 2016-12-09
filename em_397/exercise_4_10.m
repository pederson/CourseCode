% Exercise 4.10
% Bayesian Inverse Problems
% D. Pederson

P = [0.9, 0.075, 0.025;
     0.15, 0.8, 0.05;
     0.25, 0.25, 0.5;];
 
 [V,D] = eig(P);
 
 bar(V(:,1)/sum(V(:,1)))