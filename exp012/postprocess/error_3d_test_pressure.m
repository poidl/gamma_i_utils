close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')

[s1,values1]=error_3d(gamma_initial,SA,CT,p);

load('../data/gamma_i.mat')
%load('../..exp010/data/plots.mat')

[s_pressure,values_pressure]=error_3d(gamma_i,SA,CT,p);


load('data/input_data.mat')
bdy= 170<=longs(:) & longs(:)<=270 & -1<=lats(:) & lats(:)<=1;

bdy_pressure=gamma_initial(bdy);

g=gamma_rf(s(bdy),ct(bdy))
g_bb=interp1(bdy_pressure,g,values_pressure);

gamma_i=f2g(gamma_i,g,22,1); % transfer to g

[s2,values2]=error_3d(gamma_i,SA,CT,p,g_bb);
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


plot(g_bb,vdiff_pressure,'b')
hold on
plot(g_bb,vdiff_pressure,'bv')
%keyboard;

save('../data/plots.mat','g_bb','vdiff_pressure')

xlim([26.2,28])
print('-dpng','-r200',['../figures/D_f_pressure.png'])

figure()
hist(gamma_i(:),50)
xlim([26.2,28])

print('-dpng','-r200',['../figures/hist.png'])


