function [sig1,sig2] =  sig_vals(SP1,pt1,p1,SP2,pt2,p2)

% sig_vals                    Potential density of two neighbouring bottles
%==========================================================================
% 	
% DESCRPTITION:
%  Computes the gamma G values of two neighbouring bottles w.r.t. the mid 
%  pressure
%
% USAGE:  
%  [gamma1,gamma2] =  sig_vals(SP1,pt1,p1,SP2,pt2,p2)
%
% INPUT :			
% SP1, SP2         bottle Practical Salinities
% pt1, pt2         bottle Potential Temperatures
% lat1, lat2       bottle lattitudes
%
% OUTPUT:		
%  sig1, sig2   bottle neutral density values
%
% AUTHOR: 
%  Guillaume Serazin                                    
%
%==========================================================================

pmid = 0.5*(p1 + p2);
sig1 = rho_from_theta(SP1,pt1,pmid) - 1000;
sig2 = rho_from_theta(SP2,pt2,pmid) - 1000;

return