close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')

%[s1,values1]=error_3d(gamma_initial,SA,CT,p);

load('../data/gamma_i.mat')


[s2,values2]=error_3d_eps(gamma_i,SA,CT,p);

vdiff2=s2;



save('../data/plots_error_eps_3d.mat','values2','vdiff2','gamma_i')

%print('-dpdf','-r200',['../figures/hist.pdf'])


