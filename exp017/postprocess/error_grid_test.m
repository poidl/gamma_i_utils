close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')
load('../data/gamma_i.mat')

err=error_grid(gamma_i,SA,CT,p);

%err(abs(err)>1e-6)=nan;
%keyboard
err(err==0)=nan;
%save_netcdf03(log10(err),'D_f','../data/D_f.nc')
save_netcdf03(log10(err),'sn2_D_f','../data/sn2_D_f.nc')
%save_netcdf03(log10(err),'sg2_D_f_','../data/sg2_D_f_.nc')
%save_netcdf03(log10(err),'sg2_D_f_mod','../data/sg2_D_f_mod.nc')
