function [SP_g,pt_g,p_g] = SP_pt_p_on_cste_gammaG(SP0,pt0,p0,SP,pt,p,basin)

%SP_pt_p_on_cste_gammaG       Compute SP pt p on a constante gammaG surface
%==========================================================================
% 
% USAGE:    
%  [SP_g,pt_g,p_g] 
%           = SP_pt_p_on_cste_gammaG(SP_i,pt_i,p_i,lat_i,SP,pt,p,lat,basin)
% 
% DESCRIPTION:
%  Compute Practical Salinity, potential temperature and pressure on a
%  a constante gamma G surface using Newton Raphson in the depth program, 
%  the starting point being any bottle of any cast and establishing the
%  reference.
% 
% INPUT:
%  SP0   =  Practical Salinity of the reference bottle          [ g kg^-1 ]
%  pt0   =  potential Temperature of the reference bottle(ITS-90) [ deg C ]
%  p0    =  sea pressure of the reference bottle                   [ dbar ]
%  lat0  =  lattitude of the reference bottle
%  SP    =  Practical Salinity                                  [ g kg^-1 ]
%  pt    =  potential Temperature (ITS-90)                        [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%  lat   =  lattitude 
%
% OUTPUT:  
%  SP    =  Practical Salinity on the gamma G surface           [ g kg^-1 ]
%  pt    =  potential Temperature on the gamma G surface  (ITS-90)[ deg C ]
%  p     =  sea pressure on the gamma G surface                    [ dbar ]
%
%AUTHOR:
% Guillaume Sérazin
%
%==========================================================================

[~,lattitude,longitude]=size(SP);

SP_g=nan(lattitude,longitude);
pt_g=nan(lattitude,longitude);
p_g=nan(lattitude,longitude);

for j=1:longitude
   
    for i=1:lattitude
                               
            %Selecting the studied cast
            SPc=SP(:,i,j);
            ptc=pt(:,i,j);
            pc=p(:,i,j);
                        
            %Finding the intersection    
            [SP_g(i,j),pt_g(i,j),p_g(i,j)] = depth_gammaG(SP0,pt0,p0,SPc,ptc,pc,basin);
     end
                
end
    






