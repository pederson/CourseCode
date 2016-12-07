function [] = draw_cov_ellipse(axeshandle, xbar, Pxx, nsigma)

% eigen decomp
[V, D] = eig(Pxx);


 % range to plot over
%------------------------------------
N = 50;
theta = 0:1/N:2*pi+1/N;

% Parametric equation of the ellipse
%----------------------------------------
state(1,:) = nsigma*D(1,1)*cos(theta); 
state(2,:) = nsigma*D(2,2)*sin(theta);

% Coordinate transform (since your ellipse is axis aligned)
%----------------------------------------
X = V*state;
X(1,:) = X(1,:) + xbar(1);
X(2,:) = X(2,:) + xbar(2);

% Plot
%----------------------------------------
axes(axeshandle)
hold on
plot(axeshandle, X(1,:),X(2,:));
