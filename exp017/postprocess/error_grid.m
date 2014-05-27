function err=error_grid(va,s,ct,p)

    load('data/dy.mat')


    user_input;
    [nz,ny,nx]=size(s);

    % interpolate dx and dy onto grid points 
    if nx~=1
        dx=0.5*(dx(:,1:end-1)+dx(:,2:end));
        dx=horzcat(dx(:,end), dx); % sloppy 

    end
    dx=repmat(permute(dx,[3 1 2]),[nz 1 1 ]);
    dy=0.5*(dy(1:end-1,:)+dy(2:end,:));
    dy=vertcat(dy(end,:), dy); % sloppy
    dy=repmat(permute(dy,[3 1 2]),[nz 1 1 ]);


    [vx,vy,vz]=gradient(va,p,dx,dy);

    [sx,sy]=error_horz(s,ct,p,dx,dy);
    
    vx=vx./vz; % slope of va: x-component
    vy=vy./vz; % slope of va: y-component
    
    sdx= vx-sx;
    sdy= vy-sy;
    
    if nx~=1
        s2=sdx.^2+sdy.^2;
    else
        s2=sdy.^2;
    end
    
    K=1e3;
    err=K*s2;
    
end

function   [sx,sy]=error_horz(s,ct,p,dx,dy)
    
    user_input;
    
    [nz,ny,nx]=size(s);
    
    %%%
    sx2=circshift(s,[0 0 -1]);
    sx1=circshift(s,[0 0 1]);
    ctx2=circshift(ct,[0 0 -1]);
    ctx1=circshift(ct,[0 0 1]);


    drho_x2 = gsw_rho(sx2,ctx2,p)-gsw_rho(s,ct,p); % centered to the right of gridpoint
    drho_x1 = gsw_rho(s,ct,p)-gsw_rho(sx1,ctx1,p); % centered to the left of gridpoint
    
    ex=0.5*(drho_x1+drho_x2)./dx;
    
    if ~zonally_periodic
        ex(:,:,1)=nan;
        ex(:,:,end)=nan;
    end

    %%%
    sy2=circshift(s,[0 -1 0]);
    sy1=circshift(s,[0 1 0]);
    cty2=circshift(ct,[0 -1 0]);
    cty1=circshift(ct,[0 1 0]);
    
    drho_y2 = gsw_rho(sy2,cty2,p)-gsw_rho(s,ct,p); % centered to the right of gridpoint
    drho_y1 = gsw_rho(s,ct,p)-gsw_rho(sy1,cty1,p); % centered to the left of gridpoint
    
    ey=0.5*(drho_y1+drho_y2)./dy;

    ey(:,1,:)=nan;
    ey(:,end,:)=nan; 
    
    rho=gsw_rho(s,ct,p);
    [n2,pmid]=gsw_Nsquared(s,ct,p);
    n2=reshape(n2,[nz-1,ny,nx]);
    n2=0.5*(n2(1:end-1,:,:)+n2(2:end,:,:));
    n2=cat(1,n2,nan*ones(1,ny,nx));
    n2=cat(1,nan*ones(1,ny,nx),n2);
    
    fac=(1/9.81)*rho.*n2;
    
    sx=ex./fac;
    sy=ey./fac;
    
end

function [vx,vy,vz]=gradient(va,p,dx,dy)

    user_input;
    
    dp=circshift(p,[-1 0 0])-p;
    dp(end,:,:)=dp(1,:,:);
    if any(dp(1)*ones(size(dp))~=dp)
        error('grid has non-equidistant spacing in z: not implemented')
    end

    % centered fd
    vx2=circshift(va,[0 0 -1]);
    vx1=circshift(va,[0 0 1]);
    vx=(vx2-vx1)./(2*dx);


    vy2=circshift(va,[0 -1 0]);
    vy1=circshift(va,[0 1 0]);
    vy=(vy2-vy1)./(2*dy);

    vz2=circshift(va,[-1 0 0]); % z and p are oriented oppositely, as are dp and dz
    vz1=circshift(va,[1 0 0]);
    vz=(vz2-vz1)./(2*dp);
    
    if ~zonally_periodic
        vx(:,:,1)=nan;
        vx(:,:,end)=nan;
    end
    vy(:,1,:)=nan;
    vy(:,end,:)=nan;
    vz(1,:,:)=nan;
    vz(end,:,:)=nan;

end