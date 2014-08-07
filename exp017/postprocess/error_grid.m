function err=error_grid(va,s,ct,p)

    load('data/dy.mat')


    user_input;
    [nz,ny,nx]=size(s);
    
    if 0
        % interpolate dx and dy onto grid points 
        if nx~=1
            dx=0.5*(circshift(dx,[0 1])+dx);
        end

        dy=0.5*(circshift(dy,[1 0])+dy);

        dx=repmat(permute(dx,[3 1 2]),[nz 1 1 ]);
        dy=repmat(permute(dy,[3 1 2]),[nz 1 1 ]);

        [vx,vy,vz]=gradient(va,p,dx,dy);
        vx=vx./vz; % slope of va: x-component
        vy=vy./vz; % slope of va: y-component

        [sx,sy]=error_horz(s,ct,p,dx,dy);

        sdx= vx-sx;
        sdy= vy-sy;     
    else
        
        [sdx,sdy]=error_iso(va,s,ct,p,dx,dy);
    end
    
    sdx2=sdx(:).^2;
    sdy2=sdy(:).^2;

    ix=~isnan(sdx2) & isnan(sdy2);
    iy=isnan(sdx2) & ~isnan(sdy2);
    both=~isnan(sdx2) & ~isnan(sdy2);
    
    s2=nan*sdx;
    s2(ix)=sdx2(ix);
    s2(iy)=sdy2(iy);
    s2(both)=sdx2(both)+sdy2(both);
    
    K=1e3;
    err=K*s2;
    
end


function [sdx,sdy]=error_iso(va,s,ct,p,dx,dy)
    
    user_input
    
    [nz,ny,nx]=size(s);
    
    sdx=nan*ones(size(s));
    sdy=nan*ones(size(s));
    
    % north 
    sn=circshift(s,[0 -1 0]);
    ctn=circshift(ct,[0 -1 0]);
    pn=circshift(p,[0 -1 0]);    
    van=circshift(va,[0 -1 0]); 
      
    % east
    se=circshift(s,[0 0 -1]);
    cte=circshift(ct,[0 0 -1]);
    pe=circshift(p,[0 0 -1]);    
    vae=circshift(va,[0 0 -1]);
    
    for kk=1:nz
    %for kk=2:2
    %keyboard
        vsurf=squeeze(va(kk,:,:));
        
        ssurf=squeeze(s(kk,:,:));
        ctsurf=squeeze(ct(kk,:,:));
        psurf=squeeze(p(kk,:,:));
        
        if nx==1
            ssurf=ssurf';
            ctsurf=ctsurf';
            psurf=psurf';
        end
        
        vpx=var_on_surf_stef(pe,vae,vsurf);
        vsx=(vpx-psurf)./dx;        
        
        vpy=var_on_surf_stef(pn,van,vsurf);
        vsy=(vpy-psurf)./dy;

        [tr,tr,npx]=depth_ntp_simple(ssurf(:)',ctsurf(:)',psurf(:)',se(:,:),cte(:,:),pe(:,:),0*ssurf(:)');
        npx=reshape(npx,[ny,nx]);
        nsx=(npx-psurf)./dx;
        
        [tr,tr,npy]=depth_ntp_simple(ssurf(:)',ctsurf(:)',psurf(:)',sn(:,:),ctn(:,:),pn(:,:),0*ssurf(:)');
        npy=reshape(npy,[ny,nx]);
        nsy=(npy-psurf)./dy;
        sdx(kk,:,:)=vsx-nsx;
        sdy(kk,:,:)=vsy-nsy;  
        %sdx(kk,:,:)=npx;
        %sdy(kk,:,:)=npy;   
        %keyboard
       
    end

    if ~zonally_periodic
        sdx(:,:,end)=nan;
    end
    sdy(:,end,:)=nan;
    
end


function [drho]=dens_diff(va,surf,s,ct,p, ss, cts, ps)

        s1=var_on_surf_stef(s,va,surf);
        ct1=var_on_surf_stef(ct,va,surf);
        p1=var_on_surf_stef(p,va,surf);
        
        pmid= 0.5*(ps+p1);

        drho=gsw_rho(s1,ct1,pmid)-gsw_rho(ss,cts,pmid);
        


end

% function [sdx,sdy]=error_iso(va,s,ct,p,dx,dy)
% 
%     [nz,ny,nx]=size(s);
%     
%     sdx=nan*ones(size(s));
%     sdy=nan*ones(size(s));
%     
%     % north 
%     sn=circshift(s,[0 0 -1]);
%     ctn=circshift(ct,[0 0 -1]);
%     pn=circshift(p,[0 0 -1]);    
%     van=circshift(va,[0 0 -1]); 
%     
%     % south
%     ss=circshift(s,[0 0 1]);
%     cts=circshift(ct,[0 0 1]);
%     ps=circshift(p,[0 0 1]);    
%     vas=circshift(va,[0 0 1]);    
% 
%     % east
%     se=circshift(s,[0 -1 0]);
%     cte=circshift(ct,[0 -1 0]);
%     pe=circshift(p,[0 -1 0]);    
%     vae=circshift(va,[0 -1 0]);
%     
%     % west
%     sw=circshift(s,[0 1 0]);
%     ctw=circshift(ct,[0 1 0]);
%     pw=circshift(p,[0 1 0]);    
%     vaw=circshift(va,[0 1 0]);
%     
%     [n2,pmid]=gsw_Nsquared(s,ct,p);
%     n2=reshape(n2,[nz-1,ny,nx]);
%     pmid=reshape(pmid,[nz-1,ny,nx]);
%     
%     for kk=1:nz
%     %for kk=20:20    
%         surf=squeeze(va(kk,:,:));
%         
%         ssurf=squeeze(s(kk,:,:));
%         ctsurf=squeeze(ct(kk,:,:));
%         psurf=squeeze(p(kk,:,:));
%         
%         drho_x1=dens_diff(van,surf,sn,ctn,pn, ssurf, ctsurf, psurf);
%         drho_x2=-dens_diff(vas,surf,ss,cts,ps, ssurf, ctsurf, psurf);
%         ex=0.5*(drho_x1+drho_x2)./dx;
%         
%         drho_y1=dens_diff(vae,surf,se,cte,pe, ssurf, ctsurf, psurf);
%         drho_y2=dens_diff(vaw,surf,sw,ctw,pw, ssurf, ctsurf, psurf);
%         ey=0.5*(drho_y1+drho_y2)./dy;
%         
%         rhosurf=gsw_rho(ssurf,ctsurf,psurf);
%         n2surf=var_on_surf_stef(n2,pmid,psurf);
%         
%         
%         fac=(1/9.81)*rhosurf.*n2surf;
%     
%         sdx(kk,:,:)=ex./fac;
%         sdy(kk,:,:)=ey./fac;
%    
%     end
% end



% function   [sx,sy]=error_horz(s,ct,p,dx,dy)
%     
%     user_input;
%     
%     [nz,ny,nx]=size(s);
%     
%     %%%
%     sx2=circshift(s,[0 0 -1]);
%     sx1=circshift(s,[0 0 1]);
%     ctx2=circshift(ct,[0 0 -1]);
%     ctx1=circshift(ct,[0 0 1]);
% 
% 
%     drho_x2 = gsw_rho(sx2,ctx2,p)-gsw_rho(s,ct,p); % centered to the right of gridpoint
%     drho_x1 = gsw_rho(s,ct,p)-gsw_rho(sx1,ctx1,p); % centered to the left of gridpoint
%     
%     ex=0.5*(drho_x1+drho_x2)./dx;
%     
%     if ~zonally_periodic
%         ex(:,:,1)=nan;
%         ex(:,:,end)=nan;
%     end
% 
%     %%%
%     sy2=circshift(s,[0 -1 0]);
%     sy1=circshift(s,[0 1 0]);
%     cty2=circshift(ct,[0 -1 0]);
%     cty1=circshift(ct,[0 1 0]);
%     
%     drho_y2 = gsw_rho(sy2,cty2,p)-gsw_rho(s,ct,p); % centered to the right of gridpoint
%     drho_y1 = gsw_rho(s,ct,p)-gsw_rho(sy1,cty1,p); % centered to the left of gridpoint
%     
%     ey=0.5*(drho_y1+drho_y2)./dy;
% 
%     ey(:,1,:)=nan;
%     ey(:,end,:)=nan; 
%     
%     rho=gsw_rho(s,ct,p);
%     [n2,pmid]=gsw_Nsquared(s,ct,p);
%     n2=reshape(n2,[nz-1,ny,nx]);
%     n2=0.5*(n2(1:end-1,:,:)+n2(2:end,:,:));
%     n2=cat(1,n2,nan*ones(1,ny,nx));
%     n2=cat(1,nan*ones(1,ny,nx),n2);
%     
%     fac=(1/9.81)*rho.*n2;
%     
%     sx=ex./fac;
%     sy=ey./fac;
%     
% end
% 
% function [vx,vy,vz]=gradient(va,p,dx,dy)
% 
%     user_input;
%     
%     dp=circshift(p,[-1 0 0])-p;
%     dp(end,:,:)=dp(1,:,:);
%     if any(dp(1)*ones(size(dp))~=dp)
%         error('grid has non-equidistant spacing in z: not implemented')
%     end
% 
%     % centered fd
%     vx2=circshift(va,[0 0 -1]);
%     vx1=circshift(va,[0 0 1]);
%     vx=(vx2-vx1)./(2*dx);
% 
% 
%     vy2=circshift(va,[0 -1 0]);
%     vy1=circshift(va,[0 1 0]);
%     vy=(vy2-vy1)./(2*dy);
% 
%     vz2=circshift(va,[-1 0 0]); % z and p are oriented oppositely, as are dp and dz
%     vz1=circshift(va,[1 0 0]);
%     vz=(vz2-vz1)./(2*dp);
%     
%     if ~zonally_periodic
%         vx(:,:,1)=nan;
%         vx(:,:,end)=nan;
%     end
%     vy(:,1,:)=nan;
%     vy(:,end,:)=nan;
%     vz(1,:,:)=nan;
%     vz(end,:,:)=nan;
% 
% end