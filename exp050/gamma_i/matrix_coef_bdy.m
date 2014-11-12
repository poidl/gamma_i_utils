function coeff=matrix_coef_bdy(bdy,n_bdy)

w_bdy=ones(n_bdy,1);

% get N2
%keyboard
%[n2,~]=n2_smooth(s,ct,p);
%n2(n2(:)<=1e-6)=1e-6;

%[n2,pmid]=gsw_Nsquared(s,ct,p);
%n2=cat(1,n2,n2(end,:,:));
%w_bdy=1./n2(bdy);

coeff=w_bdy.*ones(n_bdy,1);