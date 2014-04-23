clear all
close all
% add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../../../omega/ansu_utils/external_scripts/'))
addpath(genpath('.'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lat=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lat');
% lon=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lon');
% p=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','pressure');
% s=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','s');
% t=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','t');
% gamma=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','gamma');
% 
% lat=lat(2:end-1); % cut what is not in gk_ak_gamma.mat
% s=s(:,:,2:end-1);
% t=t(:,:,2:end-1);
% gamma=gamma(:,:,2:end-1);
% 
% [nz nx ny]=size(s);
% lat=repmat(lat,[1 nx]); lat=repmat(lat,[1 1 nz ]);
% lon=repmat(lon,[1 ny]); lon=repmat(lon,[1 1 nz ]);
% lat=permute(lat,[3 1 2]);
% lon=permute(lon,[3 2 1]);
% 
% p=repmat(p,[1,nx]); p=repmat(p,[1 1 ny]);
% 
% s=permute(s, [1 3 2]);
% t=permute(t, [1 3 2]);
% p=permute(p, [1 3 2]);
% gamma_96=permute(gamma, [1 3 2]);
% 
% s=gsw_SA_from_SP(s,p,lon,lat);
% ct=gsw_CT_from_t(s,t,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('data/input_data.mat')
load('data/input_data_idealized.mat')
% this is to remove nan points floating in the ocean:
kinds=1:size(s,1); 
for ii=1:size(s,3)
    for jj=1:size(s,2)
        igood=~isnan(s(:,jj,ii));
        if any(igood)
            kg=find(igood,1,'first'); % kg may not be equal to zero: sea ice
            deeper=(kinds' > kg);
            kk=find(isnan(s(:,jj,ii)) & deeper,1,'first'); % shallowest nan below data
            s(kk:end,jj,ii)=nan;
            ct(kk:end,jj,ii)=nan;
%        p(kk:end,jj,ii)=nan;
        end
    end
end
lon=longs;
lat=lats;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_i = gamma_3d(s,ct,p,lon,lat);






