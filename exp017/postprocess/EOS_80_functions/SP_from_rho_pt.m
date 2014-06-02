function SP = SP_from_rho_pt(rho,pt,p)

% gsw_SP_from_rho_pt                         Absolute Salinity from density
% =========================================================================
%
% USAGE:
%  SP = gsw_SP_from_rho_pt(rho,pt,p), or equivalently

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
   error('gsw_SP_from_rho_pt:  Requires three inputs')
end %if

[md,nd] = size(rho);
[mt,nt] = size(pt);
[mp,np] = size(p);

if (mt ~= md | nt ~= nd)
    error('gsw_SP_from_rho_pt: rho and pt must have same dimensions')
end

if (mp == 1) & (np == 1)               % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)          % p is row vector,
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (np == 1)          % p is column vector,
    p = p(:,ones(1,nd));               % copy across each row.
elseif (nd == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                            % transposed then
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_SP_from_rho_pt: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho = rho.';
    pt = pt.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------



v_lab = ones(size(rho))./rho;
rho_0=rho_from_theta(zeros(size(rho)),pt,p);
rho_50=rho_from_theta(50*ones(size(rho)),pt,p);
v_0 =  ones(size(rho))./rho_0;
v_50 =  ones(size(rho))./rho_50;
 
SP = 50*(v_lab - v_0)./(v_50 - v_0);            % initial estimate of SA.

[Ior] = find(SP < 0 | SP > 50);
if ~isempty(Ior)
  SP(Ior) = NaN;
end

v_SP = (v_50 - v_0)./50; %initial estimate of v_SP, the SP derivative of v

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------

for Number_of_iterations = 1:2
    
    SP_old = SP;
    delta_v =  ones(size(rho))./rho_from_theta(SP_old,pt,p) - v_lab;
    SP = SP_old - delta_v./v_SP ; % this is half way through the modified N-R method
    SP_mean = 0.5*(SP + SP_old);
    [rho,rho_s,~,~]=eosall_from_theta(SP_mean,pt,p);
    beta=1./rho.*rho_s;
    v_SP = - beta./rho; 
    SP = SP_old - delta_v./v_SP;
    [Ior] = find(SP < 0 | SP > 50);
    if ~isempty(Ior)
        SP(Ior) = NaN; 
    end
end

% After two iterations of this modified Newton-Raphson iteration,
% the error in SP is no larger than 8x10^-13 g/kg, which 
% is machine precision for this calculation. 
 


 
if transposed
    SP = SP.';
end

end