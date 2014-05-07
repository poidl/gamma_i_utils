function err = error_3d(va,sa,ct,p)

user_input;

[nz,ny,nx]=size(va);
% va=va(:,:);
% sa=sa(:,:);
% ct=ct(:,:);
% p=p(:,:);

% find deepest location
pn=p; pn(isnan(sa))=nan;
imax=max(pn(:));
ibdy=int16(imax)/nz;

bdy=va(ibdy,:);
bdy_top=bdy(1);
bdy_bottom=bdy(end);

ps=var_on_surf_stef(p,va,bdy(20)*ones(ny,nx));
ss=var_on_surf_stef(sa,p,ps);
cts=var_on_surf_stef(ct,p,ps);

[ex,ey] = delta_tilde_rho(ss,cts,ps); % ex and ey are defined on staggered grids
% regrid
ex=0.5*(ex+circshift(ex,[0 1]));
if ~zonally_periodic % could extrapolate?
    ex(:,1)=nan;
    ex(:,end)=nan;
end
ey=0.5*(ey+circshift(ey,[1 0]));
ey(1,:)=nan;
ey(end,:)=nan;


vap=ex.^2+ey.^2;
h=imagesc(vap);
set(h,'alphadata',~isnan(vap));
set(gca,'YDir','normal');
colorbar()

figure()
vap=ex
h=imagesc(vap);
set(h,'alphadata',~isnan(vap));
set(gca,'YDir','normal');
colorbar()

figure()
vap=ey
h=imagesc(vap);
set(h,'alphadata',~isnan(vap));
set(gca,'YDir','normal');
colorbar()
end



