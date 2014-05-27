function [SPns,ptns,pns,iter] = e_solve(SP,pt,p,e,k,SP0,pt0,p0)

% e_solve                                     find the zero of the function
%==========================================================================
% 
% USAGE:         [SPns,ptns,pns,iter] = e_solve(SP,pt,p,e,k,SP0,pt0,p0)
%
% DESCRIPTION:
%  Find the zero of the e function using a bisection method, designed for
%  neutral tangent plane
%   
% INPUT:
%  SP	=	array of cast Practical Salinities                     [ g/kg ]
%  pt	=	array of cast Potential Temperatures                  [ deg C ]
%  p	=	array of cast pressures                                  [ db ]
%  e	=	array of cast e values
%  k	=	interval (k-1, k) contains the zero  
%  SP0  =	the bottle Practical Salinity                          [ g/kg ]
%  pt0  =	the bottle Potential Temperature                      [ deg C ]
%  p0   =	the bottle pressure                                      [ db ]
%
% OUTPUT:	
%  SPns  =  Practical Salinity of e zero                           [ g/kg ]
%  ptns  =	Potential Temperature of e zero                       [ deg C ]
%  pns	 =  pressure of e zero                                       [ db ]       
%
% AUTHOR: 
%  DRJ on 17/06/03
%  Revised by Guillaume Sérazin 
%  
%==========================================================================


kl=k(1);
ku=k(2);
pl = p(kl);
el = e(1);
pu = p(ku);
eu = e(2);

iter = 0;
success = 0;

while success == 0
    iter = iter + 1;
    pm = 0.5*(pl + pu);
    [SPm, ptm] = stp_interp([SP(kl),SP(ku)],[pt(kl),pt(ku)],[p(kl),p(ku)],pm);
    [sigl, sigu] = sig_vals(SP0,pt0,p0,SPm,ptm,pm);
    em = sigu - sigl;
    if el*em < 0
        pu = pm;
        eu = em;
    elseif em*eu < 0
        pl = pm;
        el = em;
    elseif em == 0
        SPns = SPm;
        ptns = ptm;
        pns = pm;
        success = 1;
    end
    if success == 0
        if (abs(em) <= 5e-5) && (abs(pu-pl) <= 5e-3 && pm>0)
            SPns = SPm;
            ptns = ptm;
            pns = pm;
            success = 1;
        elseif iter <= 10
            success = 0;
        else
            disp('WARNING in e-solve')
            disp(['iter: ', int2str(iter), '  em: ', num2str(em), '  dp: ', num2str(abs(pu-pl))])
            SPns = NaN;
            ptns = NaN;
            pns = NaN;
            success = 0;
            return
        end
    end
end


return