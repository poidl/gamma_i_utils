close all
clear all

%df = ncread('../exp017/data/D_f.nc','D_f');
%dfg = ncread('../sandbox/Guillaume/data/Df_guillaume_2d.nc','Df_guillaume_2d');
df = ncread('../exp017/data/sn2_D_f.nc','sn2_D_f');
%df = ncread('../sandbox/Guillaume/data/sn2_Df_guillaume_2d_.nc','sn2_Df_guillaume_2d_');
dfg = ncread('../sandbox/Guillaume/data/sn2_Df_guillaume_2d.nc','sn2_Df_guillaume_2d');
%df = ncread('../exp017/data/sg2_D_f.nc','sg2_D_f');
%dfg = ncread('../sandbox/Guillaume/data/sg2_Df_guillaume_2d.nc','sg2_Df_guillaume_2d');
%dfg = ncread('../exp017/data/sg2_D_f_mod.nc','sg2_D_f_mod');

df=permute(df,[3 2 1]);
dfg=permute(dfg,[3 2 1]);

v1=df;
v2=dfg;

setnan=isnan(df(:))|isnan(dfg(:));
v1(setnan)=nan;
v2(setnan)=nan;

cmin=min([df(:);dfg(:)]);
cmax=max([df(:);dfg(:)]);
%cmin=min(dfg(:));
%cmax=max(dfg(:));
%cmin=min(df(:));
%cmax=max(df(:));

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])
va=squeeze(v1);
h=imagesc(va);
set(h,'alphadata',~isnan(va))
%set(gca,'YDir','normal')
colorbar()
caxis([cmin cmax])
title('df')
print('-dpdf','-r400',['figures/df2d.pdf'])

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])
va=squeeze(v2);
h=imagesc(va);
set(h,'alphadata',~isnan(va))
%set(gca,'YDir','normal')
colorbar()
caxis([cmin cmax])
title('dfg')
print('-dpdf','-r400',['figures/dfg2d.pdf'])

figure()
va=squeeze(v1-v2);
h=imagesc(va);
set(h,'alphadata',~isnan(va))
%set(gca,'YDir','normal')
colorbar()
