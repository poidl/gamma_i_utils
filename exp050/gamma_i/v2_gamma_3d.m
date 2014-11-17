function gamma_i = v2_gamma_3d(s,ct,p,lon,lat)

% Written by D.R. Jackett
% Modified by P.M. Barker (2014)
% Modified by S. Riha (2014)
% Principal Investigators: T.J. McDougall, D.R. Jackett

[nz,ny,nx] = size(s);

user_input;

tic
if 1
    make_intersections(s,ct,p);      
end
load('./data/intersections.mat');
display(['Runtime spent on root finding: ',num2str(toc),' seconds'])

wet=~isnan(s); %

[irow,jcol,n_lateral,j_e_l,j_w_l,j_n_l,j_s_l]=matrix_ij_lateral(wet,k_east,k_west,k_north,k_south);

ibb=backbone_index(squeeze(lon(1,:,:)),squeeze(lat(1,:,:)));
[ibdy,jbdy,n_bdy,bdy]=matrix_ij_bdy(wet,ibb);

jcol=[jcol;jbdy];
irow=[irow;n_lateral+ibdy];

coeff=matrix_coef_lateral(wet,k_east,k_west,k_north,k_south,...
                              r_east,r_west,r_north,r_south,...
                              j_e_l, j_w_l, j_n_l,  j_s_l,...
                              n_lateral);
                      
coeff_bdy=matrix_coef_bdy(bdy,n_bdy);
coeff=[coeff;coeff_bdy];

A = sparse(irow,jcol,coeff,n_lateral+n_bdy,sum(wet(:)));
gamma_initial=v2_initial_guess(s,ct,p,lon,lat,ibb);
y_bdy=rhs_bdy(gamma_initial,bdy,wet,n_bdy);
y=[zeros(n_lateral,1); y_bdy];

gamma_initial=gamma_initial(wet);

if 1
    disp('starting LSQR()')
    tic
    [gamma,flag,relres,iter,resvec,lsvec] = lsqr(A,y,1e-15,10000,[],[],gamma_initial(:));
    display(['LSQR() took ',num2str(toc),' seconds for ',num2str(length(lsvec)),' iterations']);

    if length(lsvec)==length(resvec)
        mynorm=lsvec./resvec;
    else
        mynorm=lsvec./resvec(2:end);
    end
    disp(['Arnorm/(anorm*rnorm) final: ', num2str(mynorm(end))])
    disp(['Flag: ', num2str(flag)])
    save('data/mynorm.mat','mynorm')
else
    gamma = (A'*A)\(A'*y);
end

cr=get_coupled_regions(A);
disp(['Number of decoupled regions (incl. backbone region): ',num2str(length(cr))])
idecoupled=[];
for ii=1:length(cr)
    reg=cr{ii};
    if ~ismember(ibb,reg) % will not work if there is sea ice
        idecoupled=[idecoupled;reg];
    end
end
disp(['Number of grid points decoupled from backbone: ', num2str(length(idecoupled))])

gamma(idecoupled)=nan; % setting decoupled points to nan

gamma_i=nan*wet;
gamma_i(wet)=gamma;
gamma_i=reshape(gamma_i,[nz,ny,nx]);

% store some useful fields

save('data/gamma_i.mat','gamma_i')
save_netcdf03(gamma_i,'gamma_i','data/gamma_i.nc')

isdecoupled=nan*zeros(nz,ny,nx); % make nc file indicating decoupled points
dec=nan*ones(size(gamma));
dec(idecoupled)=1;
isdecoupled(wet)=dec;
save_netcdf03(isdecoupled,'isdecoupled','data/isdecoupled.nc')

save_netcdf03(s,'s','data/s.nc')

vars={'gamma_i','gamma_initial','A','y'};
save('data/lsqr_input.mat',vars{:});

end





