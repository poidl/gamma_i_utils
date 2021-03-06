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
    load('../../exp016/data/gamma_initial.mat')
    load('../../exp016/data/gamma_i.mat')
    s=SA;
    ct=CT;
    gamma=gamma_i;
end

badinds=[74,76];
s(:,:,badinds)=nan;
ct(:,:,badinds)=nan;
gamma(:,:,badinds)=nan;

% inds=(73:90);
% s=s(:,:,inds);
% ct=ct(:,:,inds);
% p=p(:,:,inds);
% gamma=gamma(:,:,inds);
% lat=lat(:,:,inds);
% lon=lon(:,:,inds);

%[Df,se_ntp,sn_ntp] = diffusivity_gamma_ref(s,ct,p,gamma,lat,lon,se_ntp,sn_ntp)
tic
[Df_guillaume,se_ntp,sn_ntp] = diffusivity_gamma_ref(s,ct,p,gamma,lat,lon);
display(['it took ',num2str(toc),' seconds']);

save_netcdf03(log10(Df_guillaume),'Df_guillaume','data/sg2_guillaume.nc')
