clear all; clear globals;
o = matrixfree;
m = o.m;
n = o.n;
r = 128;

Z = rangeA(@(x)o.A(x), @(x)o.At(x), m, n, r);
