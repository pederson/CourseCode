function [y_xt] = wave_1D_bounded(x, t, alpha, f_0, N)
% x: vector of x locations for which to solve the wave equation
% t: vector of time steps
% alpha: wave speed
% f_0: initial waveform -- must have same dimensions as x
% N: number of modes to include in the expansion
%
% y_xt: matrix with each row containing a time step of the solution
% 
% Solves by summation of Cn*sin(k*x)cos(alpha*k*t)


kn = (1:N)*pi/x(end);
dx = x(2)-x(1);


% calculate the coefficients Cn by numerical quadrature
Cn = zeros(size(kn));
for k=1:length(kn)
    num=0;
    denom=0;
    for i=1:length(x)
       num = num + f_0(i)*sin(kn(k)*x(i));
       denom = denom + (sin(kn(k)*x(i))^2);
    end
    Cn(k) = num/denom;
end

% sum up the series
y_xt = zeros(length(t), length(x));
for ti=1:length(t)
    for i=1:length(x)
        for n=1:length(Cn)
            y_xt(ti, i) = y_xt(ti, i) + Cn(n)*sin(kn(n)*x(i))*cos(alpha*kn(n)*(ti));
        end
    end
end

    