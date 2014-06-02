close all
clear all

addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('../../../david_woolsey'))
load('../exp017/postprocess/hoit.mat')


pmid=0.5*(s0+s);

s0_=s0*ones(size(s));
ct0_=ct0*ones(size(s));

F=gsw_rho(s0_,ct0_,pmid)-gsw_rho(s,ct,pmid);


plot(F,-p)

disp('Method             SA       CT       P')
[SAns,CTns,pns] = depth_ntp_jackett(s0,ct0,p0,s',ct',p');
disp(['Jackett:           ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_guillaume(s0,ct0,p0,s',ct',p');
disp(['Guillaume:         ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_jackett_fzero(s0,ct0,p0,s',ct',p');
disp(['depth_ntp_iter:    ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])
[SAns,CTns,pns] = depth_ntp_simple(s0,ct0,p0,s,ct,p,0*s0);
disp(['depth_ntp_simple:  ',num2str(SAns), '  ', num2str(CTns), '  ', num2str(pns)])

