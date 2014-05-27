function [SP_g,pt_g,p_g] = SP_pt_p_on_cste_gamma_ref(SP0,pt0,p0,gamma0,SP,pt,p,gamma)

%SP_pt_p_on_cste_gamma_ref    Compute SP pt p on a constante gamma variable
%==========================================================================
% 
% USAGE:    
%  [SP_g,pt_g,p_g] 
%           = SP_pt_p_on_cste_gamma_ref(SP0,pt0,p0,gamma0,SP,pt,p,gamma)
% 
% DESCRIPTION:
%  Compute Practical Salinity, potential temperature and pressure on a
%  a constante sigma 2 surface using Newton Raphson in the depth program, 
%  the starting point being any bottle of any cast and establishing the
%  reference.
% 
% INPUT:
%  SP0   =  Practical Salinity of the reference bottle          [ g kg^-1 ]
%  pt0   =  potential Temperature of the reference bottle(ITS-90) [ deg C ]
%  p0    =  sea pressure of the reference bottle                   [ dbar ]  
%  SP    =  Practical Salinity                                  [ g kg^-1 ]
%  pt    =  potential Temperature (ITS-90)                        [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%
% OUTPUT:  
%  SP_g  =  Practical Salinity on the gamma surface             [ g kg^-1 ]
%  pt_g  =  potential Temperature on the gamma surface(ITS-90)    [ deg C ]
%  p_g   =  sea pressure on the gamma surface                      [ dbar ]
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
   
    parfor i=1:lattitude
                               
            %Selecting the studied cast
            SPc=SP(:,i,j);
            ptc=pt(:,i,j);
            pc=p(:,i,j);
            gammac=gamma(:,i,j);
            SPc_m=min(SPc);
            SPc_M=max(SPc);
            ptc_m=min(ptc);
            ptc_M=max(ptc);
    
            %Finding the intersection          
            [SP_gamma,pt_gamma,p_gamma] = depth_gamma_ref(SP0,pt0,p0,gamma0,SPc,ptc,pc,gammac);
            
            if (p_gamma>0&&p_gamma<6500&&pt_gamma<ptc_M&&pt_gamma>ptc_m&&SP_gamma<SPc_M&&SP_gamma>SPc_m)
                SP_g(i,j)=SP_gamma;
                pt_g(i,j)=pt_gamma;
                p_g(i,j)=p_gamma;
            else
                SP_g(i,j)=NaN;
                pt_g(i,j)=NaN;
                p_g(i,j)=NaN;
            end
            
     end
                
end
    






