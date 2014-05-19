close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')

[s1,values1]=error_3d(gamma_initial,SA,CT,p);

load('../data/gamma_i.mat')

[s_pressure,values_pressure]=error_3d(gamma_i,SA,CT,p);


[s2,values2]=error_3d(gamma_i,SA,CT,p);
K=1e3;
vdiff1=K*s1;
vdiff2=K*s2;
vdiff_pressure=K*s_pressure;

% plot(values1,vdiff1,'k')
% hold on
% plot(values1,vdiff1,'ko')
plot(values2,vdiff2,'r')
hold on
plot(values2,vdiff2,'ro')
title('D_f [m^s/s]')
xlabel('\gamma^{rf} (black), \gamma^{i} (red)')


% plot(g_bb,vdiff_pressure,'b')
% hold on
% plot(g_bb,vdiff_pressure,'bv')
% 

xlim([26.2,28])
print('-dpng','-r200',['../figures/D_f.png'])

figure()
hist(gamma_i(:),50)
xlim([26.2,28])

%print('-dpdf','-r200',['../figures/hist.pdf'])


