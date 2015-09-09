function div_vel = divergence(U,V)
% Dylan Pederson
% This function calculates the divergence field given the velocity 
% components U and V assuming that the grid goes from 0 to 2pi in 
% both directions and the fields are periodic. 

% preliminaries
N = length(U);
H = 2*pi/N;
j=1:N;
jm1 = j-1;
jm1(jm1<1) = N-1;
jp1 = j+1;
jp1(jp1>N) = 2;

% approximate the derivative of u in the x direction
Du_x = (U(:,jp1) - U(:,jm1))/(2*H);
Dv_y = (V(jp1,:) - V(jm1,:))/(2*H);

div_vel = Du_x + Dv_y;
end