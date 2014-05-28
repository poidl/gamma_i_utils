close all
clear all

%df = ncread('../exp016/data/sg2_D_f.nc','D_f');
%dfg = ncread('../sandbox/Guillaume/data/sg2_guillaume.nc','Df_guillaume');
df = ncread('../exp017/data/D_f.nc','D_f');
dfg = ncread('../sandbox/Guillaume/data/Df_guillaume_2d.nc','Df_guillaume_2d');

df=permute(df,[3 2 1]);
dfg=permute(dfg,[3 2 1]);

v1=df;
v2=dfg;

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