close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')
s=SA;
ct=CT;
%load('../data/gamma_i.mat')

load('/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp624/data/gamma_i_ribs.mat')

[nz,ny,nx]=size(s);
ns=size(pns3d,1);

lon1=squeeze(lon(1,1,:));
lat1=squeeze(lat(1,:,1));
blon=187<=lon1 & lon1<=189; % 16 S, 188 E
blat=-17<=lat1 & lat1<=-15;
if (sum(blon)~=1 || sum(blat)~=1)
    error('something is wrong')
end
ilon=find(blon);
ilat=find(blat);

rpns3d=nan*pns3d;
p_s=nan*pns3d;
rp_bb=nan*ones(ns,1);
p_bb=nan*ones(ns,1);

for kk=1:ns
%    gins=var_on_surf_stef(gamma_i,p,pns3d(kk,:,:));
%     gi_bb=gins(ilat,ilon);
%     gins=gins-gi_bb;
%     gins3d(kk,:,:)=gins;

    p_bb(kk)=pns3d(kk,ilat,ilon);
    
    ptmp=p_bb(kk)+0*s;
    rp=gsw_rho(s,ct,ptmp);
    rpns3d(kk,:,:)=var_on_surf_stef(rp,p,pns3d(kk,:,:)); 
    
    rp_bb(kk)=rpns3d(kk,ilat,ilon);
    
    
    rp_s=rp_bb(kk)*ones(ny,nx);
    p_s(kk,:,:)=var_on_surf_stef(p,rp,rp_s);
    % for branched iso-surfaces (vertical inversions) this may have gone
    % wrong. check:
    if abs(p_s(kk,ilat,ilon)-pns3d(kk,ilat,ilon))>1e-9
        %error('stop')
        keyboard
    end
end

save('../data/rpns3d.mat','rpns3d','pns3d','p_s','rp_bb','p_bb','ilat','ilon')
save_netcdf03(rpns3d,'rpns3d','../data/rpns3d.nc')


dp_rp_3d=nan*pns3d;
for kk=1:ns
    if kk==93
        keyboard
   dp_rp_3d(kk,:,:)= p_s(kk,:,:)-pns3d(kk,:,:);
    end
end  
save_netcdf03(dp_rp_3d,'dp_rp_3d','../data/dp3_rp_d.nc')
%save_netcdf03(p_s,'p_s','../data/p_s.nc')

