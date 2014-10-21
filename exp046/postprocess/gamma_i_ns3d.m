close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')
s=SA;
ct=CT;
load('../data/gamma_i.mat')

load('/home/nfs/z3439823/mymatlab/omega/ansu_utils/exp631/data/gamma_i_ribs.mat')

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

gins3d=nan*pns3d;
p_s=nan*pns3d;
gi_bb=nan*ones(ns,1);
p_bb=nan*ones(ns,1);

for kk=1:ns
%    gins=var_on_surf_stef(gamma_i,p,pns3d(kk,:,:));
%     gi_bb=gins(ilat,ilon);
%     gins=gins-gi_bb;
%     gins3d(kk,:,:)=gins;
    gins3d(kk,:,:)=var_on_surf_stef(gamma_i,p,pns3d(kk,:,:)); 
    
    gi_bb(kk)=gins3d(kk,ilat,ilon);
    p_bb(kk)=pns3d(kk,ilat,ilon);
    
    gi_s=gi_bb(kk)*ones(ny,nx);
    p_s(kk,:,:)=var_on_surf_stef(p,gamma_i,gi_s);
    % for branched iso-surfaces (vertical inversions) this may have gone
    % wrong. check:
    if abs(p_s(kk,ilat,ilon)-pns3d(kk,ilat,ilon))>1e-10
        error('stop')
        %keyboard
    end
end

save('../data/gins3d_s2_no_n2.mat','gins3d','pns3d','p_s','gi_bb','p_bb','ilat','ilon')
save_netcdf03(gins3d,'gins3d','../data/gins3d.nc')


dp3d=nan*pns3d;
for kk=1:ns
   dp3d(kk,:,:)= p_s(kk,:,:)-pns3d(kk,:,:);
end  
save_netcdf03(dp3d,'dp3d','../data/dp3d.nc')
save_netcdf03(p_s,'p_s','../data/p_s.nc')

