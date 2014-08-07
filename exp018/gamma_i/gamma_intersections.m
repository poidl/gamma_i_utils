function [k_out,r_out] = gamma_intersections(s,ct,p,shift)

% Written by D.R. Jackett
% Modified by P.M. Barker (2014)
% Modified by S. Riha (2014)
% Principal Investigator: T.J. McDougall

[nz,ny,nx]=size(s);
nxy=nx*ny;

s=s(:,:);
ct=ct(:,:);
p=p(:,:);

s_=circshift(s,[0,shift]);
ct_=circshift(ct,[0,shift]);
p_=circshift(p,[0,shift]);

k_out=nan*ones(nz,nxy);
r_out=nan*ones(nz,nxy);

for ii=1:nxy
    if mod(ii,500)==0
        disp(['ii= ',num2str(ii)]);
    end
    
    if any(~isnan(s(:,ii))) && any(~isnan(s_(:,ii)))

        s_e=repmat(s_(:,ii),[1,nz]);
        ct_e=repmat(ct_(:,ii),[1,nz]);
        p_e=repmat(p_(:,ii),[1,nz]);

        [dum1,dum2,p_n]=depth_ntp_simple(s(:,ii)',ct(:,ii)',p(:,ii)', s_e, ct_e, p_e, 0*p_e);

        inan=isnan(p_n);
        p_n_stacked=repmat(p_n,[nz 1]);
        k_e=sum(p_n_stacked>=p_e,1);
        ibot= k_e==nz;
        k_e(inan|ibot)=1; % dummy
        k_e3d= k_e+nz*[0:nz-1];
        p1=p_e(k_e3d);
        p2=p_e(k_e3d+1);
        
        r=(p_n-p1)./(p2-p1);
        
        r(ibot)=0;    % Lower cast bottle is undefined!
        k_e(ibot)=nz; % Lower cast bottle is undefined!

        r(inan)=nan;
        k_e(inan)=nan;

        k_out(:,ii)=k_e;
        r_out(:,ii)=r;
        
    else
        k_out(:,ii)=nan;
        r_out(:,ii)=nan;
    end    
    
end

k_out=reshape(k_out,[nz,ny,nx]);
r_out=reshape(r_out,[nz,ny,nx]);




