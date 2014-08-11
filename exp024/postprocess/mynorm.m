close all
clear all

load('data/mynorm.mat');
semilogy(mynorm)
xlabel(iterations)
ylabel('$\frac{a}{b}$','interpreter','latex')

