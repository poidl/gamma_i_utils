clear all
close all
% add library paths
restoredefaultpath
addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../../../omega/ansu_utils/external_scripts/'))
addpath(genpath('../../../omega/ansu_utils/exp644/')) % for backbone_index
%addpath(genpath('/home/z3439823/mymatlab/gamma')) % for calling gamma_n()
addpath(genpath('../../../gamma_n'))
addpath(genpath('.'))

tic
[s,ct,p,gamma_96,lats,longs]=gammanc_to_sctp;

[s,ct,gamma_96]=make_mask_96(s,ct,p,gamma_96,longs,lats);
longs=double(longs);
lats=double(lats);
vars = {'s','ct','p','lats','longs'};
save('data/input_data.mat',vars{:})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if 1
%     sinusoidal
% else
%     subsample_gk_ak_gamma
% 
% 
%     load('data/input_data_idealized.mat')
% 
%     remove nan points floating in the ocean:
%     kinds=1:size(s,1); 
%     for ii=1:size(s,3)
%         for jj=1:size(s,2)
%             igood=~isnan(s(:,jj,ii));
%             if any(igood)
%                 kg=find(igood,1,'first'); % kg may not be equal to zero: sea ice
%                 deeper=(kinds' > kg);
%                 kk=find(isnan(s(:,jj,ii)) & deeper,1,'first'); % shallowest nan below data
%                 s(kk:end,jj,ii)=nan;
%                 ct(kk:end,jj,ii)=nan;
%            p(kk:end,jj,ii)=nan;
%             end
%         end
%     end
% 
%     vars = {'s','ct','p','lats','longs'};
%     save('data/input_data.mat',vars{:})
% end
%zero_hel_no_continents_width1

% the lats/longs are only used to calculate epsilon in error_3d()
[nz,ny,nx]=size(s);

la=squeeze(lats(1,:,:));
lo=squeeze(longs(1,:,:));

if nx==1
    la=la';
    lo=lo';
end

[dy,dx]=scale_fac(la,lo);
save('data/dy.mat', 'dx','dy') 

% % remove all except largest region
% regions=find_regions(squeeze(s(1,:,:)));
% %keyboard
% [zi,yi,xi]=size(s);
% imaxregion=0; % index of largest region
% npmaxreg=0; % number of points in largest region
% for i=1:length(regions)
%     if length(regions{i})>npmaxreg;
%         npmaxreg=length(regions{i});
%         imaxreg=i;
%     end
% end
% setnan=true(1,xi*yi);
% setnan(regions{imaxreg})=false;
% 
% s(:,setnan)=nan;
% ct(:,setnan)=nan;

va=squeeze(ct(:,:,1));
%h=imagesc(lats(1,:,1),p(:,1,1),va);
h=imagesc(va);
set(h,'alphadata',~isnan(va))

% create barrier
%s(6:end,15,:)=nan;
%ct(6:end,15,:)=nan;


% fill isolated depressions
%s(77:end,37,:)=nan;
%ct(77:end,37,:)=nan;
%nans=true(size(s));
%nans(:,22:23,:)=false;
%nans(4:end,:,:)=true;
% nans(1:8+50,:,:)=true;
% nans(13+50:end,:,:)=true;
% nans(1:8+60,:,:)=true;
% nans(12+60:end,:,:)=true;
% nans(12+59,:,:)=true; % create an isolated depression
%s(nans)=nan;
%ct(nans)=nan;

plt_transect

%keyboard

lon=longs;
lat=lats;

gamma_i = gamma_3d(s,ct,p,lon,lat);

display(['Total runtime ',num2str(toc),' seconds'])
