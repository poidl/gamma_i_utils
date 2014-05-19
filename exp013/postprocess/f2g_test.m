close all
clear all

addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('.'))

load('data/input_data.mat')
bdy= 170<=longs(:) & longs(:)<=270 & -1<=lats(:) & lats(:)<=1;
g=gamma_rf(s(bdy),ct(bdy))

load('data/gamma_i.mat')

gamma_new=f2g(gamma_i,g,22,1);


figure()


pp=squeeze(p(:,1,1));
lat=squeeze(lats(1,:,1));

va=gamma_new;
vals=plt_get_values(va(:),5,10);
vals_label=num2str(vals.','%.2f');
vals=str2num(vals_label);

va=squeeze(va(:,:,1));
%h=imagesc(lats(1,:,1),pp(:,1,1),va);
h=imagesc(lat,pp,va);
set(h,'alphadata',~isnan(va))
hold on 
[c,h]=contour(lat,pp,va,vals,'color','k');
clabel(c,h,'labelspacing',1e10)

colorbar()
title('?')

%print('-dpng','-r200',['figures/test.png'])

keyboard