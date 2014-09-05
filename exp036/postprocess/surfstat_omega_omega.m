% plotting script
clear all;
close all;

load('/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp631/data/gamma_i_ribs.mat')
va=pns3d;

%'gins3d','pns3d','p_s','gi_bb','p_bb','ilat','ilon'
load('../data/gins3d.mat') 



%compare the height of two omega surfaces

[ns,ny,nx]=size(va);

if 0
    for kk=1:ns
        va(kk,:,:)=va(kk,:,:)-va(kk,ilat,ilon);
    end
else
    [ns,ny,nx]=size(va);
    for kk=1:ns
        va(kk,:,:)=va(kk,:,:)-pns3d(kk,:,:);
        if abs(va(kk,ilat,ilon))>1e-6
            keyboard
            error('problem')
        end
    end
end

y=nan*ones(ns,1);
for kk=1:ns
    y(kk)=max(abs(va(kk,:)));
%     if y(kk)>1600
%         keyboard
%     end
end

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

plot(p_bb,y)
hold on
plot(p_bb,y,'o','markersize',3)
grid on
%ylabel('maximal variation of \gamma^i on \omega')
%ylabel('maximum of p_{\gamma^i}-p_{\omega} [dbar]')
%ylabel('maximal variation of \sigma_{p_{backbone}} on \omega [kg/m^3]')
ylabel('maximum of p_{\omega,s^2}-p_{\omega,d\rho} [dbar]')
xlabel('pressure at backbone [dbar]')

print('-dpdf','-r200',['../figures/surfstat_dp_omeg_omeg.pdf'])

