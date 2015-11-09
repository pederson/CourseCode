function [A, b, epsilon, xstar] = problem( i)

load data;

A = data(i).A;
b = data(i).b;
epsilon = data(i).pert;
xstar = data(i).xstar;
