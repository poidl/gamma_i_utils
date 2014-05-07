% create ocean with zero helicity
clear all
close all
infile='data/input_data.mat';
load(infile);
xi_new=3;
[zi,yi,xi]=size(s);

ii=find(longs(1,1,:)>337,1,'first');

% va=squeeze(ct(:,:,ii));
% h=imagesc(lats(1,:,1),p(:,1,1),va);
% set(h,'alphadata',~isnan(va))

ct2=repmat(ct(:,:,ii),[1 1 xi_new]);
s2=repmat(s(:,:,ii),[1 1 xi_new]);

p=p(:,:,1:xi_new);
lats=lats(:,:,1:xi_new);

longs=[0:178:356];
%longs=(270);
longs=repmat(longs,[yi 1]);
longs=repmat(permute(longs,[3,1,2]),[zi 1 1]);


ct=ct2;
s=s2;

outfile=infile;
vars = {'s','ct','p','lats','longs'};
save(outfile,vars{:})
