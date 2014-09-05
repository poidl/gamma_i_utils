function [c1,ce1,ce2,cw1,cw2,cn1,cn2,cs1,cs2]=weighting_coeffs_N2(n2,gam,j1,i_e,i_e_lower,i_w,i_w_lower,i_n,i_n_lower,i_s,i_s_lower)


n2=n2(gam);
%n2(isnan(n2))=nanmean(n2); % n2 is not necessarily well defined where gamma is. quick fix for testing.
if any(~isfinite(n2(:)))
    %keyboard
    error('problem')
end

c1=(1./n2(j1));
%keyboard
ce1=c1(i_e);
ce2=c1(i_e_lower);

cw1=c1(i_w);
cw2=c1(i_w_lower);

cn1=c1(i_n);
cn2=c1(i_n_lower);

cs1=c1(i_s);
cs2=c1(i_s_lower);

end