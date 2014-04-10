function gamma_i = gamma_3d(SA,CT,p,lon,lat)


% Written by D.R. Jackett
% Modified by P.M. Barker (2014)
% Modified by S. Riha (2014)
% Principal Investigator: T.J. McDougall


%Initial estimate of the neutral surface - it would be better to use a
%locally referenced density surface.
gamma_initial = gamma_rf(SA,CT); % (SA,CT,p,lon,lat);

% save_netcdf03(gamma_initial,'gamma_initial','gamma_initial.nc')
% save_netcdf03(gamma_initial-gamma_96,'gamma_diff','gamma_diff.nc')

[nz,ny,nx] = size(SA);

zonally_periodic=true;

% %load gk_interp_gamma_boundary
% [I_bg, gamma_bdry] = gamma_boundary_gammas(gamma_initial,lon,lat);

write=false;
if 1
    % east
    [k_east,r_east] = gamma_intersections(SA,CT,p,-ny);
    if ~zonally_periodic
        k_east(:,:,end)=nan;
        r_east(:,:,end)=nan;
    end
    if write
        vars = {'k_east', 'r_east'};
        save('./data/intersections_east.mat', vars{:});
    end

    % west
    [k_west,r_west] = gamma_intersections(SA,CT,p,ny);
    if ~zonally_periodic
        k_west(:,:,1)=nan;
        r_west(:,:,1)=nan;
    end
    if write        
        vars = {'k_west', 'r_west'};
        save('./data/intersections_west.mat', vars{:});
    end
    
    % north
    [k_north,r_north] = gamma_intersections(SA,CT,p,-1);
%    keyboard
    k_north(:,end,:)=nan;
    r_north(:,end,:)=nan;
    if write    
        vars = {'k_north', 'r_north'};
        save('./data/intersections_north.mat', vars{:});
    end
    
    % south
    [k_south,r_south] = gamma_intersections(SA,CT,p,1);
    k_south(:,1,:)=nan;
    r_south(:,1,:)=nan;
    if write   
        vars = {'k_south', 'r_south'};
        save('./data/intersections_south.mat', vars{:});   
    end
else    
    load('./data/intersections_east.mat');
    load('./data/intersections_west.mat');
    load('./data/intersections_north.mat');
    load('./data/intersections_south.mat');    
end

gam=~isnan(gamma_initial(:));
east=~isnan(k_east(:));
west=~isnan(k_west(:));
north=~isnan(k_north(:));
south=~isnan(k_south(:));

%notiso =  east | west | north | south; % true where gamma is not isolated; at least one eqn.

% numbering well definied gammas
sreg=cumsum(gam);

[j_e,j_e_lower]=get_jcols(sreg,-nz*ny,k_east, east);
[j_w,j_w_lower]=get_jcols(sreg, nz*ny,k_west, west);
[j_n,j_n_lower]=get_jcols(sreg,   -nz,k_north,north);
[j_s,j_s_lower]=get_jcols(sreg,    nz,k_south,south);

gamma_initial=gamma_initial(gam);
east = east(gam);
west = west(gam);
north = north(gam); 
south = south(gam);

nox=sum(gam); % number of unknowns
no_eq= east+west+north+south; % number of lateral equations at grid point 
neq=sum(no_eq); % number of lateral equations

i1=1:neq; % row index of matrix coef. 1
j1=nan*ones(neq,1); % column index of matrix coef. 1

i_e=nan*ones(1,length(j_e)); % row index of matrix coef. -(1-r)
i_w=nan*ones(1,length(j_w));
i_n=nan*ones(1,length(j_n)); 
i_s=nan*ones(1,length(j_s)); 

jstart=1;
c_e=1; % counter east
c_w=1; % counter west
c_n=1; % counter north
c_s=1; % counter south

for ix=1:nox

    jend=jstart+no_eq(ix)-1;
    j1(jstart:jend)=ix;
    cnt=0;
    if east(ix)
        i_e(c_e)=jstart+cnt;
        cnt=cnt+1;
        c_e=c_e+1;
    end
    if west(ix)
        i_w(c_w)=jstart+cnt;
        cnt=cnt+1;   
        c_w=c_w+1;
    end
    if north(ix)
        i_n(c_n)=jstart+cnt;
        cnt=cnt+1;  
        c_n=c_n+1;
    end
    if south(ix)
        i_s(c_s)=jstart+cnt;
        cnt=cnt+1; 
        c_s=c_s+1;        
    end
    jstart=jend+1;    
end
%gam
i_e_lower=i_e; % row index of matrix coef. -r
i_w_lower=i_w;
i_n_lower=i_n;
i_s_lower=i_s;


% boundary
bdy= 170<=lon(:) & lon(:)<=270 & -10<=lat(:) & lat(:)<=10;
bdy= gam & bdy;
j1_bdy= sreg(bdy); % column indices for matrix coef. 1
i1_bdy=(neq+1:neq+sum(bdy));
neq_lateral=neq;
neq_total=neq+sum(bdy);
    

irow=[i1,i_e,i_e_lower,...
         i_w,i_w_lower,...
         i_n,i_n_lower,...
         i_s,i_s_lower,... 
         i1_bdy];
  
jcol=[j1;j_e;j_e_lower;...
         j_w;j_w_lower;...
         j_n;j_n_lower;...
         j_s;j_s_lower;...
         j1_bdy];
     
 
coeff=[ones(neq_lateral,1); -(1-r_east(east));  r_east(east); ...
                         -(1-r_west(west));  r_west(west); ...    
                         -(1-r_north(north)); r_north(north); ...
                         -(1-r_south(south)); r_south(south); ...
        ones(sum(bdy),1)];
                         
     
A = sparse(irow,jcol,coeff,neq_total,nox);

b=[zeros(neq_lateral,1); gamma_initial(bdy)];

gamma = lsqr(A,b,1e-15,10000,[],[],gamma_initial(:));


end


function [j,j_lower]=get_jcols(sreg,shift,k,good)

    sreg_shifted=circshift(sreg,shift);
    k_flat=flatten3d(k);
    k_flat(~good)=1; % dummy
    j=sreg_shifted(k_flat); % column index of matrix coef. -(1-r)
    j=j(good);
    if min(j(:))==0
        keyboard
    end
    j_lower= j+1; % column index of matrix coef. -r

end

function var=flatten3d(var)
    var=var(:,:);
    [nz,nxy]=size(var);
    tt=nz*(0:nxy-1);
    tt=repmat(tt,[nz,1]);  
    var=var+tt;
    var=var(:);  
end

