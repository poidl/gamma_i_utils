function g3d=f2g_new(f,g,j,i)

% f: 3d scalar field
% i: zonal index of backbone cast
% j: meridional index of backbone cast
% g: value of g at backbone cast

[nz,ny,nx] = size(f);

fmin=min(f(:));
fmax=max(f(:));
fbb=f(:,j,i); % cast of f at backbone

% extend
df_ext_bot=fbb(end-1)-fbb(end);
dg_ext_bot=g(end-1)-g(end);
while fbb(end)<fmax
    fbb=[fbb; fbb(end)-df_ext_bot];
    g=[g; g(end)-dg_ext_bot];
end
df_ext_surf=fbb(1)-fbb(2);
dg_ext_surf=g(1)-g(2);

while fbb(1)>fmin
    fbb=[fbb(1)+df_ext_surf; fbb];
    g=[g(1)+dg_ext_surf; g];
end
nzbb=length(fbb);

g3d=nan*f;
for ii=1:ny*nx
    g3d(:,ii)=interp1(fbb,g,f(:,ii),'pchip');
end

end
