close all
clear all
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../exp024'))
addpath(genpath('.'))

load('data/input_data_gammanc.mat')

load('../exp024/data/plots_error_3d.mat')
x2=values2;
y2=vdiff2;

[nz,ny,nx]=size(s);
s(s<-90)=nan;
ct(ct<-90)=nan;
gamma(gamma<-90)=nan;

la=squeeze(lats(1,:,:));
lo=squeeze(longs(1,:,:));

if nx==1
    la=la';
    lo=lo';
end

[dy,dx]=scale_fac(la,lo);
save('data/dy.mat', 'dx','dy') 

[s,values]=error_3d(gamma,s,ct,p);

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

K=1e3;
vdiff=K*s;

h1=plot(values,vdiff,'r')
hold on
plot(values,vdiff,'ro')
h2=plot(x2,y2,'k')
plot(x2,y2,'k*')

xl1=21;
xl2=28.5;
ylim([0 3e-6]);

xl1=26;
xl2=28.5;
ylim([0 0.15e-6])

xlim([xl1,xl2]);
ylabel('D_f [m^2/s]')
xlabel('\gamma_n')
grid on
%xlabel('\gamma^{rf} (black), \gamma^{i} (red)')
ax1=gca;
pos=get(gca,'position');
set(gca,'color','none');

histax=axes('position',pos);
xhi=linspace(xl1,xl2,100);
[n1,x1]=histc(gamma(:),xhi);
[n2,x2]=histc(gamma_i(:),xhi);
%h3=plot(xhi,n)
pp1=patch([xhi xl2 xl1],[n1' 0 0], 0.85*[1 0 0],'edgecolor','none')
hold on
pp2=patch([xhi xl2 xl1],[n2' 0 0], 0.3*[1 1 1],'edgecolor','none')
alpha(0.2);

set(histax,'xticklabel',[],'xtick',[])
set(histax,'yticklabel',[],'ytick',[])

uistack(histax,'bottom')
xlim([xl1,xl2]);

%axis off

% figure()
% hist(vdiff(:))

legend([h1 h2  pp1 pp2],'location','northwest','\gamma_{n}','\gamma_i (backbone: \gamma_{n})',...
    'frequency distribution \gamma_{n}','frequency distribution \gamma_{i}')
%legend([h1 h2 ],'backbone: \gamma_{rf}','backbone: pressure')
print('-dpng','-r200',['figures/D_f_3d_compare_local.pdf'])



