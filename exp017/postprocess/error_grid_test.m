close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')
load('../data/gamma_i.mat')

err=error_grid(gamma_i,SA,CT,p);

%err(abs(err)>1e-6)=nan;

save_netcdf03(log10(err),'D_f','../data/D_f.nc')

keyboard
