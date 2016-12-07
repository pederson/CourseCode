function [h, H] = p3_measurement(x)

H = [0, 0.5; 1, 0.5];
h = H*x;