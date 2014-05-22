close all
clear all
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../exp015'))
addpath(genpath('.'))

load('data/input_data.mat')

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



xl1=25
xl2=28;

xlim([xl1,xl2]);
ylim([0 0.15e-6])
%ylim([0 0.3e-7])
ylabel('D_f [m^2/s]')
xlabel('\gamma_n')
grid on
%xlabel('\gamma^{rf} (black), \gamma^{i} (red)')
ax1=gca;
pos=get(gca,'position');
set(gca,'color','none');

histax=axes('position',pos);
xhi=linspace(xl1,xl2,100);
[n,x]=histc(gamma(:),xhi);
%h3=plot(xhi,n)
pp=patch([xhi xl2 xl1],[n' 0 0], 0.85*[1 1 1],'edgecolor','none')

set(histax,'xticklabel',[],'xtick',[])
set(histax,'yticklabel',[],'ytick',[])

uistack(histax,'bottom')
xlim([xl1,xl2]);

%axis off

figure()
hist(vdiff(:))

legend([h1  pp],'\gamma_{n}','frequency distribution')
%legend([h1 h2 ],'backbone: \gamma_{rf}','backbone: pressure')
print('-dpdf','-r200',['figures/D_f_gamma_n.pdf'])



