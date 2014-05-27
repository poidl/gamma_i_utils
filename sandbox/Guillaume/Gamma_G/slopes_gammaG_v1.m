function [se,sn] = slopes_gammaG_v1(SP,pt,p,lat,long)

% slopes_gammaG                      Calculate the slope of Gamma G surface
%==========================================================================
% 
% USAGE:    
%  [se,sn] = slope_gammaG(SP,pt,p,lat,long)
% 
% DESCRIPTION:
%  Calculate the slopes based on the difference of pressure between the
%  intersection of a bottle in a neighbouring cast and the equivalent 
%  bottle of the cast using the program depth_gammaG. The difference of 
%  pressure is converted in a difference in height to find the slope using
%  a geometrical definition
% 
% INPUT:
%  SP    =  Practical Salinity                                  [ g kg^-1 ]
%  pt    =  potential Temperature (ITS-90)                        [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%  lat   =  lattitude 
%  long  =  longitude
%
% OUTPUT:  
%  ss  =  slope error
%  sx  =  x-component of slope 
%  sy  =  y-component of slope 
%
%AUTHOR:
% Guillaume Sérazin
%
%==========================================================================


[depth,lattitude,longitude]=size(SP);
se=nan(size(SP));
sn=se;
p_N=nan(size(SP));
p_E=nan(size(SP));
deltap_N=nan(size(SP));
deltap_E=nan(size(SP));

% Beginning of the calculation

for i=1:longitude-1
   
    disp(['gamma G: ',num2str(i)])
    
    for j=1:(lattitude-1)
        
        %Selecting the East test cast
        SPt_E=SP(:,j,i+1);
        ptt_E=pt(:,j,i+1);
        pt_E=p(:,j,i+1);
        latt_E=lat(1,j,i+1);
        longt_E=long(1,j,i+1);

        %Selecting the North test cast
        SPt_N=SP(:,j+1,i);
        ptt_N=pt(:,j+1,i);
        pt_N=p(:,j+1,i);
        latt_N=lat(1,j+1,i);
        longt_N=long(1,j+1,i);

        %Selecting the studied cast
        SP0=SP(:,j,i);
        pt0=pt(:,j,i);
        p0=p(:,j,i);
        lat0=lat(1,j,i);
        long0=long(1,j,i);
        
        parfor k=1:depth
            
            %Selecting a bottle of the studied cast
            SP0b=SP0(k);
            pt0b=pt0(k);
            p0b=p0(k);
            z0b=gsw_z_from_p(p0b,lat0);

            %Computing Easterm slope
            [~,~,p_gamma_E] = depth_gammaG_v1(SP0b,pt0b,p0b,lat0,long0,SPt_E,ptt_E,pt_E,latt_E,longt_E);
            z_gamma_E=gsw_z_from_p(p_gamma_E,latt_E);
            deltap_E(k,j,i)=p0b-p_gamma_E;
            deltaz_E=z0b-z_gamma_E;
            deltax=gsw_distance([long0 longt_E],[lat0 latt_E] ,[p0b p_gamma_E]);
            p_E(k,j,i)=p_gamma_E;
%             if(isreal(deltaz_E))
              se(k,j,i)=deltaz_E/deltax;
%             end
            

            %Computing Northern slope
            [~,~,p_gamma_N] = depth_gammaG_v1(SP0b,pt0b,p0b,lat0,long0,SPt_N,ptt_N,pt_N,latt_N,longt_N);
            z_gamma_N=gsw_z_from_p(p_gamma_N,latt_N);
            deltaz_N=z0b-z_gamma_N;
            deltap_N(k,j,i)=p0b-p_gamma_N;
            deltay=gsw_distance([long0 longt_N],[lat0 latt_N] ,[p0b p_gamma_N]);
            p_N(k,j,i)=p_gamma_N;
%             if(isreal(deltaz_N))
              sn(k,j,i)=deltaz_N/deltay;
%             end

        end
        
    end
    
end

return



   