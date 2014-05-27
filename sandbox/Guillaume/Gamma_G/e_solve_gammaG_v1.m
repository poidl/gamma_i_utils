function [SPns,ptns,pns,iter] = e_solve_gammaG_v1(SP,pt,p,lat,long,e,k,SP0,pt0,p0,lat0,long0)

% e_solve_gammaG                              find the zero of the function
%==========================================================================
% 
% USAGE:         [SPns,ptns,pns,iter] = e_solve(SP,pt,p,e,k,SP0,pt0,p0)
%
% DESCRIPTION:
%  Find the zero of the e function using a bisection method, designed for
%  gamma G
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
%  Adapted by Guillaume Sérazin for use on Gamma G
%  
%==========================================================================


global Iprob_i Iprob_j Iprob_k
kl=k(1);
ku=k(2);
pl = p(kl);
el = e(1);
pu = p(ku);
eu = e(2);
lat1=lat(1);
long1=long(1);

iter = 0;
success = 0;

while success == 0
    iter = iter + 1;
    pm = 0.5*(pl + pu);
    [SPm, ptm] = stp_interp([SP(kl),SP(ku)],[pt(kl),pt(ku)],[p(kl),p(ku)],pm);
    [gammal, gammau] = gammaG_vals_v1(SP0,pt0,p0,lat0,long0,SPm,ptm,pm,lat1,long1);
    em = gammau - gammal;
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
        break
    end
    if success == 0
        if (abs(em) <= 5e-5) && (abs(pu-pl) <= 5e-3)
            SPns = SPm;
            ptns = ptm;
            pns = pm;
            success = 1;
        elseif iter <= 10
            success = 0;
        else
            disp('WARNING in e-solve')
            disp(['iter: ', int2str(iter), '  em: ', num2str(em), '  dp: ', num2str(abs(pu-pl)),'  lat: ',num2str(lat0)])
            disp(['i: ',num2str(Iprob_i),'  j: ', num2str(Iprob_j),'  k: ', num2str(Iprob_k)])
            SPns = NaN;
            ptns = NaN;
            pns = NaN;
            return
        end
    end
end


return