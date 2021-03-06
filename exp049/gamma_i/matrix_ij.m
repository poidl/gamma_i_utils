function [irow,jcol,n_total,n_lateral,m,j_e_l,j_w_l,j_n_l,j_s_l]=matrix_ij(wet,k_east,k_west,k_north,k_south,ibb)

[nz,ny,nx]=size(k_east);

% numbering wet data points
sreg=cumsum(wet);
sreg(~wet)=nan;

% get the wet index of vertically adjacent point pair between which we do the linear
% interpolation. j_e is the index of the upper bottle.
% note that sum(east)==length(j_e)
[j_e,j_e_lower,j_e_l]=get_jcols(sreg,-nz*ny,k_east); 
[j_w,j_w_lower,j_w_l]=get_jcols(sreg, nz*ny,k_west);
[j_n,j_n_lower,j_n_l]=get_jcols(sreg,   -nz,k_north);
[j_s,j_s_lower,j_s_l]=get_jcols(sreg,    nz,k_south);

east=~isnan(k_east(:)); % we can form an eastern equation here
west=~isnan(k_west(:));
north=~isnan(k_north(:));
south=~isnan(k_south(:));

m=sum(wet(:)); % number of unknowns
no_eq= east+west+north+south; % number of lateral equations at grid point 
n_lateral=sum(no_eq); % number of lateral equations


i1=1:n_lateral; % row index of matrix coef. 1

j1=[]; % column index of matrix coef. 1
for ii=1:length(no_eq)
    j1=[j1;ii*ones(no_eq(ii),1)];
end
nu=cumsum([east,west,north,south],2); % numbering the equations at each gridpoint starting from 1
cs=cumsum(no_eq); 
nu(2:end,:)=nu(2:end,:)+bsxfun(@times,cs(1:end-1),ones(1,4));  % shift that number according to the equations of previous grid points
i_e=nu(east,1)'; % row index of matrix coef. -(1-r)
i_w=nu(west,2)';
i_n=nu(north,3)';
i_s=nu(south,4)';

i_e_lower=i_e(j_e_l); % row index of matrix coef. -r
i_w_lower=i_w(j_w_l);
i_n_lower=i_n(j_n_l);
i_s_lower=i_s(j_s_l);


% boundary
bdy=false(nz,ny,nx);
bdy(:,ibb)=true;
%keyboard
bdy= wet & bdy(:);
n_total=n_lateral+sum(bdy);

j1_bdy= sreg(bdy); % column indices for matrix coef. 1
i1_bdy=(n_lateral+1:n_lateral+sum(bdy));
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
     