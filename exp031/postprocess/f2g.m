function g3d=f2g(f,g,j,i)

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

fbb_=repmat(fbb,[1 nz*ny*nx]);

for ii=1:nx*ny
    c=f(:)';
    c_=repmat(c,[nzbb 1]);
    up=fbb_<=c_;
    kup=sum(up,1); % number of grid points above 
    
    inan=isnan(c);
   
    %outcrop= kup==0;    
    %bottom= kup==nz;
    
    kup(inan)=1; % dummy

    g_1=g(kup);
    g_2=g(kup+1);
    fbb1=fbb(kup);
    fbb2=fbb(kup+1);

    df=(c'-fbb1)./(fbb2-fbb1);
    g3d=g_1+(g_2-g_1).*df;
    
    g3d(inan)=nan;
    
end
%keyboard
g3d=reshape(g3d,[nz ny nx]);


end
