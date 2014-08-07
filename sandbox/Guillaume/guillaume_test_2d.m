close all
clear all
restoredefaultpath
addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('.'))

if 0
    load('../data/input_data_gammanc.mat')
    s=double(s);
    ct=double(ct);
    p=double(p);
    gamma=double(gamma);
    lat=lats;
    lon=longs;
else
    load('../../exp017/data/gamma_initial.mat')
    load('../../exp017/data/gamma_i.mat')
    s=SA;
    ct=CT;
    gamma=gamma_i;
end

% badinds=[74,76];
% s(:,:,badinds)=nan;
% ct(:,:,badinds)=nan;
% gamma(:,:,badinds)=nan;

% inds=(73:90);
% s=s(:,:,inds);
% ct=ct(:,:,inds);
% p=p(:,:,inds);
% gamma=gamma(:,:,inds);
% lat=lat(:,:,inds);
% lon=lon(:,:,inds);

%[Df,se_ntp,sn_ntp] = diffusivity_gamma_ref(s,ct,p,gamma,lat,lon,se_ntp,sn_ntp)
[Df_guillaume,se_ntp,sn_ntp] = diffusivity_gamma_ref_2d(s,ct,p,gamma,lat,lon);

%save_netcdf03(log10(Df_guillaume),'sg2_Df_guillaume_2d','data/sg2_Df_guillaume_2d.nc')
%save_netcdf03(log10(Df_guillaume),'sn2_Df_guillaume_2d','data/sn2_Df_guillaume_2d.nc')
save_netcdf03(log10(Df_guillaume),'Df_guillaume_2d','data/Df_guillaume_2d.nc')
