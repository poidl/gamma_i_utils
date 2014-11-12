function    coeff=matrix_coef_lateral(r_e,  r_w,  r_n,  r_s,...
                                      j_e_l,j_w_l,j_n_l,j_s_l,...
                                      n_lateral)
                          
c1=1; ce1=1; ce2=1; cw1=1; cw2=1; cn1=1; cn2=1; cs1=1; cs2=1;
%n_bdy=n_total-n_lateral;

%w_bdy=ones(n_bdy,1);
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
                              -cs1.*(1-r_s); -cs2.*r_s(j_s_l)];
      
                          