% plotting script
clear all;
close all;

%'gins3d','pns3d','p_s','gi_bb','p_bb','ilat','ilon'
load('../data/gins3d_s2_no_n2.mat')
%load('../data/rpns3d.mat')
load('../data/gamma_initial.mat')

va=gins3d;
%va=p_s;
%va=rpns3d;

[ns,ny,nx]=size(va);

if 1
    for kk=1:ns
        va(kk,:,:)=va(kk,:,:)-va(kk,ilat,ilon);
    end
else
    [ns,ny,nx]=size(va);
    for kk=1:ns
        va(kk,:,:)=va(kk,:,:)-pns3d(kk,:,:);
    end
end

y=nan*ones(ns,1);
for kk=1:ns
    y(kk)=max(abs(va(kk,:)));
end

sz=1.5*[13 10];
figure('PaperSize',sz,'PaperPosition',[0 0 sz(1) sz(2)]) 

plot(p_bb,y)
hold on
plot(p_bb,y,'o','markersize',3)
grid on
ylabel('maximal variation of \gamma^i on \omega')
%ylabel('maximum of p_{\gamma^i}-p_{\omega} [dbar]')
%ylabel('maximal variation of \sigma_{p_{backbone}} on \omega [kg/m^3]')
%ylabel('maximum of p_{\sigma_{p_{backbone}}}-p_{\omega} [dbar]')
xlabel('pressure at backbone [dbar]')

 print('-dpdf','-r200',['../figures/surfstat_gam_on_omeg_s2_no_n2.pdf'])
% print('-dpdf','-r200',['../figures/surfstat_dp_gam_omeg_s2_no_n2.pdf'])
% print('-dpdf','-r200',['../figures/surfstat_rpot_on_omeg.pdf'])
% print('-dpdf','-r200',['../figures/surfstat_dp_rpot_omeg.pdf'])

