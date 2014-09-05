clear all
close all
% add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))

[s,ct,p]=gammanc_to_sctp;
[n2,pmid]=gsw_Nsquared(s,ct,p);
%keyboard
[nz,ny,nx]=size(s);
n2=reshape(n2,[nz-1,ny,nx]);
keyboard
n2(n2<=0)=1e-9;
save_netcdf03(log10(n2),'log n2','n2_gamma96.nc')
