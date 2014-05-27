close all
clear all
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../exp015'))
addpath(genpath('.'))

load('data/input_data_gammanc.mat')


err=error_grid(gamma,s,ct,p);

%err(abs(err)>1e-6)=nan;

save_netcdf03(log10(err),'D_f','data/D_f_gammanc.nc')

keyboard
