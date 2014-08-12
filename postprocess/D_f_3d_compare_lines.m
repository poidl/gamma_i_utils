close all
clear all
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../exp024'))
addpath(genpath('.'))


runs=[26,27,28,31,29,30];
runs=[runs,24];
cnt=0;
for ii=1:length(runs)
    cnt=cnt+1;
    fname=['../exp0',num2str(runs(ii)),'/data/plots_error_3d.mat'];
    load(fname);
    x{ii}=values2;
    y{ii}=vdiff2;
end


sz=1.5*[16 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

hh=nan*(1:length(runs));

cnt=1;
col=get(gcf,'DefaultAxesColorOrder');
for ii=1:length(runs)
    hh(ii)=semilogy(x{ii},y{ii},'color',col(cnt,:))
    if cnt==size(col,1)
        cnt=1
    else
        cnt=cnt+1;
    end      
    hold on
    %plot(x{ii},y{ii},'ro')
end


xl1=26;
xl2=28.5;
% ylim([0 3e-6]);
xlim([xl1,xl2]);

ylabel('D_f [m^2/s]')
xlabel('value of iso-surface')
grid on
%xlabel('\gamma^{rf} (black), \gamma^{i} (red)')

label={'50 \gamma^i LSQR() iterations','100','200','300','400','1000','8000+'};
legend(hh,'location','northwest',label)

print('-dpdf','-r200',['figures/D_f_3d_compare_lines.pdf'])



