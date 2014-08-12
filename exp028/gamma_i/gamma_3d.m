function gamma_i = gamma_3d(SA,CT,p,lon,lat)


% Written by D.R. Jackett
% Modified by P.M. Barker (2014)
% Modified by S. Riha (2014)
% Principal Investigator: T.J. McDougall

[nz,ny,nx] = size(SA);

%Initial estimate of the neutral surface - it would be better to use a
%locally referenced density surface.
gamma_initial = gamma_rf(SA,CT); % (SA,CT,p,lon,lat);
%gamma_initial(:)=0;

%gamma_initial=p/max(p(:));

gi=gamma_initial;

save('data/gamma_initial.mat','gamma_initial','SA','CT','p','lon','lat')
%keyboard
%sum(~isnan(gamma_initial(:)))

% l1=(1:nz)';
% l2=repmat(l1,[1,ny]);
% gamma_initial=repmat(permute(l2,[1 2 3]),[1 1 nx]);
% gamma_initial(isnan(SA))=nan;

% % analytic initial condition:
% gamma_initial=ones(nz,ny,nx);
% mid=repmat([1 2 3]',[1 3]);
% south=mid-0.4;
% north=mid+0.4;
% gamma_initial(:,1,:)=south;
% gamma_initial(:,2,:)=mid;
% gamma_initial(:,3,:)=north;
% gamma_initial(isnan(SA(:)))=nan;
% gi=gamma_initial; % this is the analytic solution;
% %gamma_initial=gi+1e-1*randn(size(gi)); % perturb analytic solution to create interesting initial condition
% gamma_initial(13:15)=gi(13:15);

%keyboard
% save_netcdf03(gamma_initial,'gamma_initial','gamma_initial.nc')
% save_netcdf03(gamma_initial-gamma_96,'gamma_diff','gamma_diff.nc')

user_input;

% %load gk_interp_gamma_boundary
% [I_bg, gamma_bdry] = gamma_boundary_gammas(gamma_initial,lon,lat);

if 0
    write=true;
    % east
    [k_east,r_east] = gamma_intersections(SA,CT,p,-ny);
    if ~zonally_periodic
        k_east(:,:,end)=nan;
        r_east(:,:,end)=nan;
    end

    % west
    [k_west,r_west] = gamma_intersections(SA,CT,p,ny);
    if ~zonally_periodic
        k_west(:,:,1)=nan;
        r_west(:,:,1)=nan;
    end
    
    % north
    [k_north,r_north] = gamma_intersections(SA,CT,p,-1);
    k_north(:,end,:)=nan;
    r_north(:,end,:)=nan;
    
    % south
    [k_south,r_south] = gamma_intersections(SA,CT,p,1);
    k_south(:,1,:)=nan;
    r_south(:,1,:)=nan;
    
    if write   
        vars = {'k_east', 'r_east','k_west', 'r_west',...
            'k_north', 'r_north','k_south', 'r_south'};
        save('./data/intersections.mat', vars{:});  
    end
    
else    
    load('./data/intersections.mat');  
end

gam=~isnan(gamma_initial(:));
east=~isnan(k_east(:));
west=~isnan(k_west(:));
north=~isnan(k_north(:));
south=~isnan(k_south(:));


%notiso =  east | west | north | south; % true where gamma is not isolated; at least one eqn.

% numbering well definied gammas
sreg=cumsum(gam);
sreg(~gam)=nan;


% [j_e,j_e_lower]=get_jcols(sreg,-nz*ny,k_east);
% [j_w,j_w_lower]=get_jcols(sreg, nz*ny,k_west);
% [j_n,j_n_lower]=get_jcols(sreg,   -nz,k_north);
% [j_s,j_s_lower]=get_jcols(sreg,    nz,k_south);
%keyboard
[j_e,j_e_lower,j_e_l]=get_jcols(sreg,-nz*ny,k_east);
[j_w,j_w_lower,j_w_l]=get_jcols(sreg, nz*ny,k_west);
[j_n,j_n_lower,j_n_l]=get_jcols(sreg,   -nz,k_north);
[j_s,j_s_lower,j_s_l]=get_jcols(sreg,    nz,k_south);


gamma_initial=gamma_initial(gam);
east = east(gam);
west = west(gam);
north = north(gam); 
south = south(gam);
r_east=r_east(gam);
r_west=r_west(gam);
r_north=r_north(gam);
r_south=r_south(gam);

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

% i_e_lower=i_e; % row index of matrix coef. -r
% i_w_lower=i_w;
% i_n_lower=i_n;
% i_s_lower=i_s;
i_e_lower=i_e(j_e_l); % row index of matrix coef. -r
i_w_lower=i_w(j_w_l);
i_n_lower=i_n(j_n_l);
i_s_lower=i_s(j_s_l);


% boundary
%bdy= 170<=lon(:) & lon(:)<=270 & -10<=lat(:) & lat(:)<=10;
%bdy= 179<=lon(:) & lon(:)<=181 & -5<=lat(:) & lat(:)<=-3; % 16 S, 188 E
bdy= 187<=lon(:) & lon(:)<=189 & -17<=lat(:) & lat(:)<=-15; % 16 S, 188 E
if sum(bdy)~=nz
    disp('WARNING: multiple backbones (or none at all)')
    keyboard
end
bdy= gam & bdy;
j1_bdy= sreg(bdy); % column indices for matrix coef. 1
i1_bdy=(neq+1:neq+sum(bdy));
neq_lateral=neq;
neq_total=neq+sum(bdy);
%keyboard

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
     
r_e=r_east(east);
r_w=r_west(west);
r_n=r_north(north);
r_s=r_south(south);
     
% r_e=r_east(east);
% r_w=r_west(west);
% r_n=r_north(north);
% r_s=r_south(south);

% coeff=[ones(neq_lateral,1); -(1-r_e); -r_e; ...
%                             -(1-r_w); -r_w; ...    
%                             -(1-r_n); -r_n; ...
%                             -(1-r_s); -r_s; ...
%         ones(sum(bdy),1)];
w_bdy=1;
coeff=[ones(neq_lateral,1); -(1-r_e); -r_e(j_e_l); ...
                            -(1-r_w); -r_w(j_w_l); ...    
                            -(1-r_n); -r_n(j_n_l); ...
                            -(1-r_s); -r_s(j_s_l); ...
          w_bdy*ones(sum(bdy),1)];
                        

A = sparse(irow,jcol,coeff,neq_total,nox);

%bdy_vals=ones(sum(bdy(gam)),1);
%b=[zeros(neq_lateral,1); bdy_vals];
load('data/gamma_96.mat')
gamma_bdy=gamma_96(gam);
b=[zeros(neq_lateral,1); w_bdy*gamma_bdy(bdy(gam))];

if 1
    %gamma = lsqr(A,b,1e-15,10000,[],[],gamma_initial(:));
    disp('starting LSQR()')
    tic
    [gamma,flag,relres,iter,resvec,lsvec] = lsqr(A,b,1e-15,200,[],[],gamma_initial(:));
    display(['LSQR() took ',num2str(toc),' seconds for ',num2str(length(lsvec)),' iterations']);
    %keyboard
    if length(lsvec)==length(resvec)
        mynorm=lsvec./resvec;
    else
        mynorm=lsvec./resvec(2:end);
    end
    disp(['Arnorm/(anorm*rnorm) final: ', num2str(mynorm(end))])
    disp(['Flag: ', num2str(flag)])
    save('data/mynorm.mat','mynorm')
else
    gamma = (A'*A)\(A'*b);
end
% set points where no equation exists to nan (their value hasn't changed from the initial condition)
su=sum(abs(A));
no_equation=su==0;
%no_equation=all(A==0);
gamma(no_equation)=nan;
gamma_i=nan*gam;
gamma_i(gam)=gamma;
gamma_i=reshape(gamma_i,[nz,ny,nx]);


save('data/gamma_i.mat','gamma_i')

gamma_i(isnan(SA(:)))=nan;
plt( squeeze(gamma_i(:,:,1)) );
% rpot=gsw_rho(SA,CT,0*p);
% plt( squeeze(rpot(:,:,1)) );

%keyboard
%cregs=find_coupled(A);


vars={'gamma_i','gamma_initial','A','b'};
save('data/lsqr_input.mat',vars{:});
keyboard
end

function cregs=find_coupled(A)
    error('this doesn''t work I think?')
    [m,n]=size(A);
    cregs=nan*ones(n,1);

    done_r=[];
    list_r=1; % row list
    ireg=1;
    keyboard
    while any(isnan(cregs))
        while ~isempty(list_r);
            r=list_r(1);
            cols=find(A(r,:)~=0);
            cregs(cols)=ireg;
            for c=cols
                rows=find(A(:,c)~=0);
                list_r=[list_r;rows];
            end   
            list_r=unique(list_r);            
            done_r=[done_r,r];
            list_r=setdiff(list_r,done_r);
        end
        keyboard
        A=A(~done_r,:);
        done_r=[];
        list_r=1; % row list
        ireg=ireg+1;
    end 
end

% function [j,j_lower]=get_jcols(sreg,shift,k)
%     good=~isnan(k(:));
%     sreg_shifted=circshift(sreg,shift);
%     k_flat=flatten3d(k);
%     k_flat(~good)=1; % dummy
%     j=sreg_shifted(k_flat); % column index of matrix coef. -(1-r)
%     j=j(good);
%     if min(j(:))==0
%         keyboard
%     end
%     j_lower= j+1; % column index of matrix coef. -r
%     outside=find(j_lower==length(good)+1);
%     if ~isempty(outside)
%         keyboard
%         j_lower(outside)=1; % dummy
%     end 
% end

function [j,j_lower,good_lower]=get_jcols(sreg,shift,k)
    [nz,ny,nx]=size(k);
    good=~isnan(k(:));
    sreg_shifted=circshift(sreg,shift);
    k_flat=flatten3d(k);
    k_flat(~good)=1; % dummy
    j=sreg_shifted(k_flat); % column index of matrix coef. -(1-r)
    j=j(good);
    
    sreg_shifted_vert=circshift(sreg_shifted,-1);
    sreg_shifted_vert(nz:nz:nz*ny*nx)=nan;

    j_lower=sreg_shifted_vert(k_flat);
    j_lower=j_lower(good);
    good_lower=~isnan(j_lower);
    
    j_lower=j_lower(good_lower);
end

function var=flatten3d(var)
    var=var(:,:);
    [nz,nxy]=size(var);
    tt=nz*(0:nxy-1);
    tt=repmat(tt,[nz,1]);  
    var=var+tt;
    var=var(:);  
end

function plt(va)
    figure()
    %h=imagesc(lats(1,:,1),p(:,1,1),va);
    h=imagesc(va);
    set(h,'alphadata',~isnan(va))
    colorbar()
end

