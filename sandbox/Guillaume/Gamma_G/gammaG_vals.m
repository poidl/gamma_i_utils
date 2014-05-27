function [gamma1,gamma2] =  gammaG_vals(SP1,pt1,p1,SP2,pt2,p2,basin)

% gammaG_vals                 Calculate Gamma G of two neighbouring bottles
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
%  gamma1, gamma2   bottle neutral density values
%
% AUTHOR: 
%  Guillaume Serazin                                    
%
%==========================================================================

gamma1 = gsw_gammaG(SP1,pt1,p1,basin);
gamma2 = gsw_gammaG(SP2,pt2,p2,basin);

return