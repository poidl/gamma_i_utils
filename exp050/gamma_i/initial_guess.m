function gamma_initial=initial_guess(s,ct,p,lon,lat,ibb)

[nz,ny,nx]=size(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition 

%gamma_initial = gamma_rf(s,ct); % (s,ct,p,lon,lat);

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