function [Df,se_ntp,sn_ntp] = diffusivity_gamma_ref(SP,pt,p,gamma,lat,long,se_ntp,sn_ntp)

% diffusivity_gammaG                                 Diffusivity of Gamma G      
%==========================================================================
% 
% USAGE:  
%  Df = diffusivity_gammaG(SP,pt,p,lat,long,basin)
%
% DESCRIPTION:
%  Calculate the fictitious diffusivity of gamma G function
% 
%
% INPUT:
%  SP       =  Practical Salinity                                  [ g/kg ]
%  pt       =  Potential temperature                              [ deg C ]
%  p        =  Pressure                                                [db]
%  gamma_n  =  neutral density                                  [ kg m^-3 ]
%  lat      =  lattitude                                                  
%  long     =  longitude
%  basin    =  acronym of the studied basin 
%
%  SP, pt, p, lat & long need to have the same dimensions.
% 
% OUTPUT:
%  Df  = Fictitious diffusivity                                [ m^2.s^-1 ] 
%
% AUTHOR: 
%  Guillaume S�razin
%
%==========================================================================



%Calculation of the slopes
if nargin==7
%   - Gamma ref    
        [se_gamma,sn_gamma] = slopes_gamma_ref(SP,pt,p,gamma,lat,long);   
    
else
%   - Gamma ref    
        [se_gamma,sn_gamma] = slopes_gamma_ref(SP,pt,p,gamma,lat,long);
        
%   - Neutral tangent plane    
        [se_ntp,sn_ntp] = slopes_ntp(SP,pt,p,lat,long);
end

        
%Lateral diffusivity        
K=1000;

%Squared difference of slopes
s2_e=(se_gamma-se_ntp).^2;
s2_n=(sn_gamma-sn_ntp).^2;

[depth,lattitude,longitude]=size(SP);
Df=NaN(depth,lattitude,longitude);


%Calculation of Df recombining the slopes
for k=1:depth
    for j=1:lattitude
        for i=1:longitude
            if(isnan(s2_e(k,j,i))==0)
                if(isnan(s2_n(k,j,i))==0)
                    Df(k,j,i)=K*(s2_n(k,j,i)+s2_e(k,j,i));
                else
                    Df(k,j,i)=K*s2_e(k,j,i);
                end
            else
                if(isnan(s2_n(k,j,i))==0)
                    Df(k,j,i)=K*s2_n(k,j,i);
                else
                    Df(k,j,i)=NaN;
                end
            end
                    
        end
    end
end


