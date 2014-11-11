function gamma_i = gamma_3d(s,ct,p,lon,lat)

% Written by D.R. Jackett
% Modified by P.M. Barker (2014)
% Modified by S. Riha (2014)
% Principal Investigator: T.J. McDougall

%error('there are decoupled systems of equations')

[nz,ny,nx] = size(s);

user_input;

% %load gk_interp_gamma_boundary
% [I_bg, gamma_bdry] = gamma_boundary_gammas(gamma_initial,lon,lat);
tic
if 0
    make_intersections(s,ct,p);
else    
    load('./data/intersections.mat');  
end
display(['Runtime spent on root finding: ',num2str(toc),' seconds'])

wet=~isnan(s(:)); %
east=~isnan(k_east(:));
west=~isnan(k_west(:));
north=~isnan(k_north(:));
south=~isnan(k_south(:));

east = east(wet);
west = west(wet);
north = north(wet); 
south = south(wet);
r_east=r_east(wet);
r_west=r_west(wet);
r_north=r_north(wet);
r_south=r_south(wet);

ibb=backbone_index(squeeze(lon(1,:,:)),squeeze(lat(1,:,:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[irow,jcol,n_total,n_lateral,m,j_e_l,j_w_l,j_n_l,j_s_l,bdy]=matrix_ij(wet,k_east,k_west,k_north,k_south,ibb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_e=r_east(east);
r_w=r_west(west);
r_n=r_north(north);
r_s=r_south(south);
     

c1=1; ce1=1; ce2=1; cw1=1; cw2=1; cn1=1; cn2=1; cs1=1; cs2=1;
n_bdy=n_total-n_lateral;
w_bdy=ones(n_bdy,1);
% get N2
%keyboard
%[n2,~]=n2_smooth(s,ct,p);
%n2(n2(:)<=1e-6)=1e-6;

%[n2,pmid]=gsw_Nsquared(s,ct,p);
%n2=cat(1,n2,n2(end,:,:));
%keyboard
%[c1,ce1,ce2,cw1,cw2,cn1,cn2,cs1,cs2]=weighting_coeffs_N2(n2,wet,...
%                   j1,i_e,i_e_lower,i_w,i_w_lower,i_n,i_n_lower,i_s,i_s_lower);

%w_bdy=1./n2(bdy);

coeff=[c1.*ones(n_lateral,1); -ce1.*(1-r_e); -ce2.*r_e(j_e_l); ...
                                -cw1.*(1-r_w); -cw2.*r_w(j_w_l); ...    
                                -cn1.*(1-r_n); -cn2.*r_n(j_n_l); ...
                                -cs1.*(1-r_s); -cs2.*r_s(j_s_l); ...
          w_bdy.*ones(n_bdy,1)];
                  
% keyboard
A = sparse(irow,jcol,coeff,n_total,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition 

%Initial estimate of the neutral surface - it would be better to use a
%locally referenced density surface.
%gamma_initial = gamma_rf(s,ct); % (s,ct,p,lon,lat);

% % analytic initial condition:
% gamma_initial=ones(nz,ny,nx);
% mid=repmat([1 2 3]',[1 3]);
% south=mid-0.4;
% north=mid+0.4;
% gamma_initial(:,1,:)=south;
% gamma_initial(:,2,:)=mid;
% gamma_initial(:,3,:)=north;
% gamma_initial(isnan(s(:)))=nan;
% gi=gamma_initial; % this is the analytic solution;
% %gamma_initial=gi+1e-1*randn(size(gi)); % perturb analytic solution to create interesting initial condition
% gamma_initial(13:15)=gi(13:15);


tis=gsw_t_from_CT(s(:,ibb),ct(:,ibb),p(:,ibb)); % in-situ
gbdy=get_gamma_n(s(:,ibb),tis,p(:,ibb),lon(:,ibb),lat(:,ibb));
%keyboard

% construct initial data set
igood=~isnan(gbdy);
if ~all(igood) % fill in some values at the bottom 
    dgam=diff(gbdy);
    dgam=dgam(igood); dgam=dgam(1:end-1);
    dgam=mean(dgam(end-5:end)); % take the mean of the 5 deepest values.
    g_deepest=gbdy(sum(igood));
    fill=g_deepest+dgam*(1:sum(~igood));
    gbdy(~igood)=fill;
end
gamma_initial=repmat(gbdy,[1,ny,nx]);
gamma_initial=gamma_initial(wet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b=[zeros(n_lateral,1); w_bdy.*gamma_initial(bdy(wet))];

if 1
    %gamma = lsqr(A,b,1e-15,10000,[],[],gamma_initial(:));
    disp('starting LSQR()')
    tic
    [gamma,flag,relres,iter,resvec,lsvec] = lsqr(A,b,1e-15,10000,[],[],gamma_initial(:));
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

% testing for decoupled systems
b_test_coupled=[zeros(n_lateral,1); w_bdy.*ones(sum(bdy(wet)),1)];
[gamma_,flag]=lsqr(A,b_test_coupled,1e-15,1000,[],[],0*gamma_initial(:));
idecoupled=gamma_==0;
disp(['Number of decoupled points: ',num2str(sum(idecoupled))])
gamma(gamma_==0)=nan;

% set points where no equation exists to nan (their value hasn't changed from the initial condition)
% not necessary/effective when testing for decoupled systems
% su=sum(abs(A));
% no_equation=su==0;
% %no_equation=all(A==0);
% gamma(no_equation)=nan;


gamma_i=nan*wet;
gamma_i(wet)=gamma;
gamma_i=reshape(gamma_i,[nz,ny,nx]);

save('data/gamma_i.mat','gamma_i')
save_netcdf03(gamma_i,'gamma_i','data/gamma_i.nc')

% some useful output

isdecoupled=zeros(nz,ny,nx);
isdecoupled(wet)=~gamma_;
isdecoupled(isdecoupled==0)=nan;
save_netcdf03(isdecoupled,'isdecoupled','data/isdecoupled.nc')

save_netcdf03(s,'s','data/s.nc')

vars={'gamma_i','gamma_initial','A','b'};
save('data/lsqr_input.mat',vars{:});

end


function [j,j_lower,good_lower]=get_jcols(sreg,shift,k)
    % find the wet indices of the vertically adjacent point pair between
    % which we do the linear interpolation
    % sreg and k have the shape of the original grid. 
    % sreg: indices of wet points
    % k: vertical index of the point above which the characteristic segment intersects 
    % j: array of indices of upper point
    % j_lower: array of indices of lower point (possibly smaller than j)
    % good_lower: logical array of size(j). true if j_lower exists.
    
    [nz,ny,nx]=size(k);
    good=~isnan(k(:)); % characteristic segments don't exist in every direction.
    sreg_shifted=circshift(sreg,shift);
    k_flat=flatten3d(k); % from vertical index to global index
    k_flat(~good)=1; % dummy. in next line k_flat is used as index. can't be nan.
    j=sreg_shifted(k_flat); % column index of matrix coef. -(1-r)
    j=j(good); % remove dummies
    
    sreg_shifted_vert=circshift(sreg_shifted,-1); % index of point below
    sreg_shifted_vert(nz:nz:nz*ny*nx)=nan; % points below bottom are nan

    j_lower=sreg_shifted_vert(k_flat);
    j_lower=j_lower(good); % discard index of lower point if there is no char. seg.
    
    good_lower=~isnan(j_lower); % the index of lower point could be nan, in case the char. seg. intersects exactly with upper point
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



