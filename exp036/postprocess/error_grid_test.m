close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/input_data.mat')
load('../data/gamma_i.mat')

tic
err=error_grid(gamma_i,s,ct,p);
display(['it took ',num2str(toc),' seconds']);
%err(abs(err)>1e-6)=nan;

D_f_grid=err;
save('../data/D_f_grid.mat','D_f_grid');
save_netcdf03(log10(err),'D_f','../data/D_f.nc')


