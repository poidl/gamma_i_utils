close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

load('../data/gamma_i.mat')

gamma=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','gamma');
lat=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lat');
lon=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lon');

inan= gamma(:)==-99;
gamma(inan)=nan;

gamma_n=permute(gamma, [1 3 2]);
gamma_n=double(gamma_n);


blon=187<=lon & lon<=189; % 16 S, 188 E
blat=-17<=lat & lat<=-15;

if (sum(blon)~=1 || sum(blat)~=1)
    error('something is wrong')
end
ilon=find(blon);
ilat=find(blat);

g_bb=gamma_n(:,ilat,ilon);

gamma_i=f2g_new(gamma_i,g_bb,ilat,ilon);

save('../data/gamma_i_f2g.mat','gamma_i')



