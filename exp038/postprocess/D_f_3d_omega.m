close all
clear all

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

load('../data/plots_error_3d.mat')
load('/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp624/data/gamma_i_ribs.mat')

x1=values2;
y1=vdiff2_bar;
y2=vdiff2_med;
%load('../exp014/data/plots.mat')
%x2=g_bb;
%y2=vdiff_pressure;
% load('../exp012/data/plots.mat')
% x3=g_bb;
% y3=vdiff_pressure;

% means
h1=semilogy(x1,y1,'b')
hold on
semilogy(x1,y1,'bo','markersize',3)
h2=semilogy(values2_omega,mdf,'r')
semilogy(values2_omega,mdf,'ro','markersize',3)

% medians
h3=semilogy(x1,y2,'--b')
semilogy(x1,y2,'bo','markersize',3)
h4=semilogy(values2_omega,df_med,'--r')
semilogy(values2_omega,df_med,'ro','markersize',3)
% h1=plot(x1,y1,'b')
% hold on
% plot(x1,y1,'bo','markersize',3)
% h2=plot(values2_omega,mdf,'r')
% plot(values2_omega,mdf,'ro','markersize',3)

%h2=plot(x2,y2,'r')
%plot(x2,y2,'ro')
% h3=plot(x3,y3,'g')
% plot(x3,y3,'go')

xl1=21;
xl2=28.5;
%ylim([1e-15 1e-7]);

%xl1=27;
%xl2=28.5;
%ylim([0 0.15e-6])

xlim([xl1,xl2]);
title('D_f [m^2/s]')
xlabel('value of iso-surface')
grid on
%xlabel('\gamma^{rf} (black), \gamma^{i} (red)')
ax1=gca;
pos=get(gca,'position');
set(gca,'color','none');

%histax=axes('position',pos,'visible','off');
histax=axes('position',pos);

xhi=linspace(xl1,xl2,100);
[n,x]=histc(gamma_i(:),xhi);
% %h3=plot(xhi,n)
 pp=patch([xhi xl2 xl1],[n' 0 0], 0.85*[1 1 1],'edgecolor','none')
 %pp=plot(xhi,n)
 %ylim([0 7000]);
% 


 xlim([xl1,xl2]);
 set(histax,'xticklabel',[],'xtick',[])
 set(histax,'yticklabel',[],'ytick',[])

% 
 uistack(histax,'bottom')
  set(histax,'visible','off');

%legend([h1 h2 pp],'backbone: \gamma_{rf}','backbone: pressure','frequency distribution')
%legend([h1  pp],'location','west','\gamma_i (backbone: \gamma_{n})','frequency distribution')
legend([h1 h3 h2 h4 pp],'location','northwest','\gamma_i (backbone: \gamma_{n})','median','\omega clamped at backbone','median','frequency distribution')
%legend([h1 h2 ],'backbone: \gamma_{rf}','backbone: pressure')
print('-dpdf','-r200',['../figures/D_f_3d_omega_global.pdf'])

%print('-dpdf','-r200',['../figures/hist.pdf'])
