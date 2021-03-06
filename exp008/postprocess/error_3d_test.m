close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_initial.mat')

e1=error_3d(gamma_initial,SA,CT,p);

load('../data/gamma_i.mat')

e2=error_3d(gamma_i,SA,CT,p);

plot(e1)
hold on
plot(e2)
plot(e2,'o')
keyboard