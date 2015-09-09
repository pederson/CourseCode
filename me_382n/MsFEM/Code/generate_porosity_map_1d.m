function porosity_map = generate_porosity_map_1d(grid)
% Dylan Pederson 
% Miles Gray
%
% This generates a fractal porosity map through the Weierstrass-Mandelbrot
% function (the real part). We set a=2.2 and x is the 1d variable from 0 to
% 1

% size of the input
[sx] = length(grid);

% number of k values
num_k = 40;
k = 1:num_k;

x = [1:sx]/sx;
a = 2.2;
porosity_map = zeros(sx, 1);
for i=1:sx
        porosity_map(i) = sum((sin(k.^a*pi*x(i)))./(pi*k.^a));
end

porosity_map = porosity_map./max(porosity_map);

%{
% visualize for testing
figure()
plot(x, porosity_map,'r')
%}