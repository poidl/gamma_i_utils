close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')
load('../data/gamma_i.mat')

[s_pressure,values_pressure]=error_3d(gamma_i,SA,CT,p);


load('data/input_data.mat')
%bdy= 170<=longs(:) & longs(:)<=270 & -1<=lats(:) & lats(:)<=1;
bdy= 187<=longs(:) & longs(:)<=189 & -17<=lats(:) & lats(:)<=-15; % 16 S, 188 E
ss=find(bdy);
[ks,ys,xs]=ind2sub(size(SA),ss(1));
bdy_pressure=gamma_initial(bdy);

g=gamma_rf(s(bdy),ct(bdy))
g_bb=interp1(bdy_pressure,g,values_pressure);

%keyboard
gamma_i=f2g(gamma_i,g,ys,xs); % transfer to g

K=1e3;
vdiff_pressure=K*s_pressure;

plot(g_bb,vdiff_pressure,'b')
hold on
plot(g_bb,vdiff_pressure,'bv')
title('D_f [m^s/s]')
xlabel('\gamma^{rf} (black), \gamma^{i} (red)')

%keyboard;

xlim([26.2,28])
print('-dpdf','-r200',['../figures/D_f_pressure.png'])

save('../data/plots.mat','g_bb','vdiff_pressure')

figure()
hist(gamma_i(:),50)
xlim([26.2,28])

print('-dpdf','-r200',['../figures/hist.png'])


