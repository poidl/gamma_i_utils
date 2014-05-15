close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')

[s1,values1]=error_3d(gamma_initial,SA,CT,p);

load('../data/gamma_i.mat')

load('data/input_data.mat')
bdy= 170<=longs(:) & longs(:)<=270 & -1<=lats(:) & lats(:)<=1;
g=gamma_rf(s(bdy),ct(bdy))
gamma_i=f2g(gamma_i,g,22,1); % transfer to g

[s2,values2]=error_3d(gamma_i,SA,CT,p);
K=1e3;
vdiff1=K*s1;
vdiff2=K*s2;

% plot(values1,vdiff1,'k')
% hold on
% plot(values1,vdiff1,'ko')
plot(values2,vdiff2,'r')
hold on
plot(values2,vdiff2,'ro')
title('D_f [m^s/s]')
xlabel('\gamma^{rf} (black), \gamma^{i} (red)')

%keyboard;


print('-dpdf','-r200',['../figures/D_f.pdf'])
