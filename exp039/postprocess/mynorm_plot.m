close all
clear all

load('data/mynorm.mat');

sz=1.0*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

mynorm=mynorm(1:300);
semilogy(mynorm)
xlabel('iterations')
ylabel('$\frac{|| A^T\cdot(b-A\cdot x)|| }{|| A||_F||b-A\cdot x||}$','interpreter','latex')
%ylabel('$\frac{|| A |}{| A |}$','interpreter','latex')

print('-dpdf','-r200',['figures/mynorm_initial.pdf'])