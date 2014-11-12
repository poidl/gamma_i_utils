function gamma_i = gamma_3d(s,ct,p,lon,lat)

% Written by D.R. Jackett
% Modified by P.M. Barker (2014)
% Modified by S. Riha (2014)
% Principal Investigator: T.J. McDougall

[nz,ny,nx] = size(s);

user_input;

tic
if 0
    make_intersections(s,ct,p);
else    
    load('./data/intersections.mat');  
end
display(['Runtime spent on root finding: ',num2str(toc),' seconds'])

wet=~isnan(s); %

[irow,jcol,n_lateral,j_e_l,j_w_l,j_n_l,j_s_l]=matrix_ij_lateral(wet,k_east,k_west,k_north,k_south);

ibb=backbone_index(squeeze(lon(1,:,:)),squeeze(lat(1,:,:)));
[ibdy,jbdy,n_bdy,bdy]=matrix_ij_bdy(wet,ibb);

jcol=[jcol;jbdy];
irow=[irow;n_lateral+ibdy];

east= ~isnan(k_east(:))  & wet(:); % estward equation exists
west= ~isnan(k_west(:))  & wet(:);
north=~isnan(k_north(:)) & wet(:);
south=~isnan(k_south(:)) & wet(:);

r_e=r_east(east);
r_w=r_west(west);
r_n=r_north(north);
r_s=r_south(south);

coeff=matrix_coef_lateral(r_e,  r_w,  r_n,  r_s,...
                          j_e_l,j_w_l,j_n_l,j_s_l,...
                          n_lateral);
                      
coeff_bdy=matrix_coef_bdy(bdy,n_bdy);
coeff=[coeff;coeff_bdy];


A = sparse(irow,jcol,coeff,n_lateral+n_bdy,sum(wet(:)));

gamma_initial=initial_guess(s,ct,p,lon,lat,ibb);

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

% testing for decoupled systems
b_test_coupled=[zeros(n_lateral,1); ones(sum(bdy(wet)),1)];
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

vars={'gamma_i','gamma_initial','A','y'};
save('data/lsqr_input.mat',vars{:});

end





