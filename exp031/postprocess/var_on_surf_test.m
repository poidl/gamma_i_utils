clear all
close all
% add library paths
load('../data/input_data.mat')

psurf1=squeeze(p(end,:,:));

ssurf1=var_on_surf_stef(s,p,psurf1);

ssurf2=squeeze(s(end,:,:));

dif=abs(ssurf1-ssurf2);

max(abs(dif(:)))

imagesc(dif)
colorbar

figure()
imagesc(ssurf1)
colorbar

figure()
imagesc(ssurf2)
colorbar