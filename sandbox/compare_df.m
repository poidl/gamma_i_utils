close all
clear all

%df = ncread('../exp016/data/sg2_D_f.nc','D_f');
%dfg = ncread('../sandbox/Guillaume/data/sg2_guillaume.nc','Df_guillaume');
df = ncread('../exp016/data/D_f.nc','D_f');
%dfg = ncread('../sandbox/Guillaume/data/Df_guillaume.nc','Df_guillaume');

df=permute(df,[3 2 1]);
%dfg=permute(dfg,[3 2 1]);

kk=75;
kk=40;

v1=df(kk,:,:);
%v2=dfg(kk,:,:);

%cmin=min([df(:);dfg(:)]);
%cmax=max([df(:);dfg(:)]);
cmin=min([df(:)]);
cmax=max([df(:)]);
sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])

va=squeeze(v1);
h=imagesc(va);
set(h,'alphadata',~isnan(va))
set(gca,'YDir','normal')
colorbar()
caxis([cmin cmax])
title('df')
print('-dpdf','-r400',['figures/df800_new.pdf'])

% sz=1.5*[13 10];
% figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)])
% va=squeeze(v2);
% h=imagesc(va);
% set(h,'alphadata',~isnan(va))
% set(gca,'YDir','normal')
% colorbar()
% caxis([cmin cmax])
% title('dfg')
% print('-dpdf','-r400',['figures/dfg1500_new.pdf'])