close all
clear all
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../exp024'))
addpath(genpath('.'))

run=24;
fname=['../exp0',num2str(run),'/data/D_f_grid.mat'];
load(fname);
fname=['../exp0',num2str(run),'/data/gamma_initial.mat'];
load(fname);

dfg1=D_f_grid;

run=25;
fname=['../exp0',num2str(run),'/data/D_f_grid.mat'];
load(fname);

dfg2=D_f_grid;

nz=size(dfg1,1);

va1=nan*ones(nz,1);
va2=nan*ones(nz,1);
for ii=1:nz
    va1(ii)=nanmean(dfg1(ii,:));
    va2(ii)=nanmean(dfg2(ii,:));
end

sz=1.5*[10 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

p=p(:,1,1);
h1=semilogx(va1,-p)
hold on
h2=semilogx(va2,-p,'r')

xlabel('global depth average of D_f [m^2/s]')
ylabel('depth [m]')
grid on

legend([h1 h2 ],'location','west','\gamma_i (backbone: \gamma_{n})','\gamma_i (backbone: p/p_{max})')


print('-dpdf','-r200',['figures/D_f_grid_z_vs_df.pdf'])

