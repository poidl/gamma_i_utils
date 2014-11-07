restoredefaultpath
addpath(genpath('../../../neutral_common'))
addpath(genpath('.'))

%path_to_gammanc='/home/nfs/z3439823/mymatlab/neutral_common/gamma_n/gamma_JackettMcDougall96';
%[s,ct,p,lon,lat]=get_input_jackett96(path_to_gammanc);
%fname='/home/nfs/z3439823/mymatlab/neutral_common/data/WOA13/woa13.nc';
%[s,ct,p,lon,lat]=get_input_woa13_4deg(fname);
fname='/home/nfs/z3439823/mymatlab/neutral_common/data/wghc/convert_to_netcdf.py/wghc.nc';
[s,ct,p,lon,lat]=get_input_wghc_4deg(fname);  
%fname='/home/nfs/z3439823/mymatlab/neutral_common/data/nemo/nemo.nc';
%[s,ct,p,lon,lat]=get_input_nemo_4deg(fname); 

[s,ct,p,lon,lat]=crop_to_gamma_n(s,ct,p,lon,lat);
[s,ct]=only_keep_largest_region(s,ct); % remove everything except largest region
    
    
field = gamma_3d_vertical(s,ct,p,lon,lat);
