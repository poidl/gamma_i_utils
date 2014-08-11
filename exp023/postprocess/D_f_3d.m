close all
clear all

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

load('../exp013/data/plots.mat')
x1=values2;
y1=vdiff2;
load('../exp014/data/plots.mat')
x2=g_bb;
y2=vdiff_pressure;
% load('../exp012/data/plots.mat')
% x3=g_bb;
% y3=vdiff_pressure;

h1=plot(x1,y1,'b')
hold on
plot(x1,y1,'bo')
h2=plot(x2,y2,'r')
plot(x2,y2,'ro')
% h3=plot(x3,y3,'g')
% plot(x3,y3,'go')

xl1=26.2
xl2=28;

xlim([xl1,xl2]);
ylim([0 0.15e-6])
%ylim([0 0.3e-7])
ylabel('D_f [m^s/s]')
%xlabel('\gamma^{rf} (black), \gamma^{i} (red)')
ax1=gca;
pos=get(gca,'position');
set(gca,'color','none');

histax=axes('position',pos);
xhi=linspace(xl1,xl2,100);
[n,x]=histc(gamma_i(:),xhi);
%h3=plot(xhi,n)
pp=patch([xhi xl2 xl1],[n' 0 0], 0.85*[1 1 1],'edgecolor','none')

set(histax,'xticklabel',[],'xtick',[])
set(histax,'yticklabel',[],'ytick',[])

uistack(histax,'bottom')
xlim([xl1,xl2]);

%axis off


legend([h1 h2 pp],'backbone: \gamma_{rf}','backbone: pressure','frequency distribution')
%legend([h1 h2 ],'backbone: \gamma_{rf}','backbone: pressure')
print('-dpdf','-r200',['figures/D_f_3d_new_.pdf'])

%print('-dpdf','-r200',['../figures/hist.pdf'])


