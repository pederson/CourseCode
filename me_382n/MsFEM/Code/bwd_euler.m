function [T, Y] = bwd_euler(M, K, q, dt, tspan, y0)
% Solves [M]dy/dt = K*y + q using the backward euler method

% force y0 to be a column vector
y0 = reshape(y0,length(y0),1);

% construct the non-changing M matrix
R = M - dt*K;

% construct the time-dependent g matrix
yi = y0;

Y = zeros(length(tspan), length(M));
T = tspan;
Y(1,:) = y0';
for i=2:length(tspan)
    g_i = M*yi + dt*q;
    
    % solve the linear system
    y_new = R\g_i;
    
    % update stuff
    Y(i,:) = y_new';
    yi = y_new;
end

end
    

