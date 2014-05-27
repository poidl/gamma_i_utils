function [SP0,pt0] = stp_interp(SP,pt,p,p0)

% stp_interp                                          Pressure interpolated
%==========================================================================
% 
% USAGE:         [SP0,pt0] = stp_interp(SP,pt,p,p0)
%  
% DESCRIPTION:
%  Linearly interpolate Practical Salinity and Potential temperature on a 
%  cast to a specified pressure
%
% INPUT:          
%  SP             cast Practical Salinities
%  pt             cast Potential Temperatures
%  p             cast pressures
%  p0            pressure level
% 
% OUTPUT:
%  s0            interpolated Practical Salinity
%  t0            interpolated Potential Temperature
%
%  DRJ on 17/06/03
%
%==========================================================================
%  Find the index of a scalar in a monotonically increasing array
n = length(p);
k = NaN;
if p(1) < p0 && p0 < p(n)
    inds_p = find(p >= p0);
    k = inds_p(1) - 1;
elseif p0 <= p(1)
    k = 1;
elseif p0 >= p(n)
    k = n - 1;
else
    SP0=NaN;
    pt0=NaN;
    return
end

r = (p0 - p(k))/(p(k+1) - p(k));
SP0 = SP(k) + r*(SP(k+1) - SP(k));
pt0 = pt(k) + r*(pt(k+1) - pt(k));
return