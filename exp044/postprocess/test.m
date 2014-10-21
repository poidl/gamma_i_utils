close all
clear all

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

load('../data/plots_error_eps_3d.mat')
load('/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp624/data/gamma_i_ribs.mat')
load('/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp624/data/gamma_i_ribs_values2.mat')

x1=values2;
y1=vdiff2;
%load('../exp014/data/plots.mat')
%x2=g_bb;
%y2=vdiff_pressure;
% load('../exp012/data/plots.mat')
% x3=g_bb;
% y3=vdiff_pressure;

h1=semilogy(x1,y1,'b')
hold on
semilogy(x1,y1,'bo','markersize',3)
h2=semilogy(values2_omega,rmsdrho,'r')
semilogy(values2_omega,rmsdrho,'ro','markersize',3)

xl1=27;
xl2=28.5;

xlim([xl1,xl2]);

grid on

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
 uistack(histax,'bottom');

ylabel(ax1,'$$\sqrt{\overline{\Delta\rho^2}}\quad \rm [kg/m^3]$$','interpreter','latex','fontsize',13)
%ylabel(histax,'$$\sqrt{\overline{\Delta\rho^2}}\quad \rm [kg/m^3]$$','interpreter','latex','fontsize',13)
xlabel('value of iso-surface')

%legend([h1 h2 pp],'backbone: \gamma_{rf}','backbone: pressure','frequency distribution')
%legend([h1  pp],'location','west','\gamma_i (backbone: \gamma_{n})','frequency distribution')
%legend([h1 h2 pp],'location','southwest','\gamma_i (backbone: \gamma_{n})','\omega clamped at backbone','frequency distribution')
%legend([h1 h2 ],'backbone: \gamma_{rf}','backbone: pressure')
print('-dpdf','-r200',['../figures/test.pdf'])

%print('-dpdf','-r200',['../figures/hist.pdf'])