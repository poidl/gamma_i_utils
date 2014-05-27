function gammaG = gsw_gammaG(SP,pt,p,basin)

% gsw_gamma_G                                               Neutral density
%==========================================================================
% 
% USAGE:  
%  rho_pt = gsw_rho_pt(SP,pt,lat)
%
% DESCRIPTION:
%  Calculates the approximate form of neutral density gamma G using
%  Practical Salinity and Potential Temperature, using the
%
% INPUT:
%  SP     =  Absolute Salinity                                     [ g/kg ]
%  pt     =  Conservative Temperature                             [ deg C ]
%  lat    =  latitude in decimal degrees north              [ -90 ... +90 ]
%  basin  =  acronym of the studied basin 
%
%  SP & pt need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SP & pt are MxN.
%
% OUTPUT:
%  gammmaG = Approximate form of Neutral density                [ kg m^-3 ]
%
% AUTHOR: 
%  Guillaume Sérazin
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_gamma_G:  Requires three inputs')
end %if

[ms,ns] = size(SP);
[mt,nt] = size(pt);


if (mt ~= ms || nt ~= ns)
    error('gsw_gamma_G: SP and pt must have same dimensions')
end



if ms == 1
    SP = SP';
    pt = pt';
   
    transposed = 1;
else
    transposed = 0;
end


%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% These few lines ensure that SP is non-negative.
[I_neg_SP] = find(SP < 0);
if ~isempty(I_neg_SP)
    SP(I_neg_SP) = 0;
end

%Normalization
SP=SP/42;
pt=pt/40;


s=0;
n=0;

%Recovering the right coefficients

switch basin(1)
    
    case {'S','N'} 
        
        switch basin
        
            case 'SPI'
                load ([pwd '\Data\Polynomials\SPI_P_coeff'])
                                    
            case 'NP'
                load ([pwd '\Data\Polynomials\NP_P_coeff'])
               
            case 'SI'
                load ([pwd '\Data\Polynomials\SI_P_coeff'])
                
                
            case 'NI'
                load ([pwd '\Data\Polynomials\NI_P_coeff'])
                
            case 'SA'
                load ([pwd '\Data\Polynomials\SA_P_coeff'])
                
            
            case 'SO'
                load ([pwd '\Data\Polynomials\SO_P_coeff'])
                s=1;
                
          
            case 'NA'
                load ([pwd '\Data\Polynomials\NA_P_coeff'])
               n=0;
                
            case 'NArc'
                load ([pwd '\Data\Polynomials\NArc_P_coeff'])              
                
 
            case 'SW'
                load ([pwd '\Data\Polynomials\SW_P_coeff'])
                                
                
            case 'SAnt'
                load ([pwd '\Data\Polynomials\SAnt_P_coeff'])
                s=0;
                                
            case 'NS_WORLD'
                load ([pwd '\Data\Polynomials\NS_WORLD_P_coeff'])
                
            otherwise
                disp([basin ' is not a valid input'])
                return 
        end       
       
    otherwise
        disp([basin ' is not a valid input'])
        return 

end
 

%Rebuilding the polynomial
% ------------------------
if(s==0)
    if(n==1)
        SP=SP*42;
        pt=pt*40;       
    end
    gammaG=Fit(1,3)*ones(ms,ns);
    for k=2:length(Fit)
        i=Fit(k,1);
        j=Fit(k,2);
        gammaG=gammaG+Fit(k,3)*(SP.^i.*pt.^j);         
    end
    

    
elseif(s==1)
    
    p_ref=700;
    pt_ref=2.5;
    c_pt=0.65;
    
    gammaG=Fit(1,3)*ones(ms,ns);
    for k=2:(length(Fit))
        i=Fit(k,1);
        j=Fit(k,2);
        gammaG=gammaG+Fit(k,3)*(SP.^i.*pt.^j);         
    end    
   
    gammaG_SO=Fit_p(1,3)*ones(ms,ns).*exp(-p/p_ref).*(1/2-1/2*tanh((40*pt-pt_ref)/c_pt));
    for k=2:(length(Fit_p))
        i=Fit_p(k,1);
        j=Fit_p(k,2);
        gammaG_SO=gammaG_SO+Fit_p(k,3)*(SP.^i.*pt.^j).*exp(-p/p_ref).*(1/2-1/2*tanh((40*pt-pt_ref)/c_pt));     
    end    
    
    gammaG=gammaG+gammaG_SO;
    
end

% gammaG_N=Fit_N(1,3)*ones(ms,ns);
% for k=2:length(Fit_N)
%     i=Fit_N(k,1);
%     j=Fit_N(k,2);
%     gammaG_N=gammaG_N+Fit_N(k,3)*(SP.^i.*pt.^j);         
% end

% gammaG_D=ones(ms,ns);
% for k=1:length(Fit_D)
%     i=Fit_D(k,1);
%     j=Fit_D(k,2);
%     gammaG_D=gammaG_D+Fit_D(k,3)*(SP.^i.*pt.^j);         
% end
% 
% gammaG=gammaG_N./gammaG_D;

    

%Denormalization
if(n==0)
    gammaG=20*gammaG-20;
end

if transposed
    gammaG= gammaG';
end

end

