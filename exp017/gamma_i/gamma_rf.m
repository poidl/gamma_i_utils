function neutral_density_rf = gamma_rf(SR,CT)

SP = gsw_SP_from_SR(SR);

%load gamma_coeffs_ct
gcoeffs =[    1.002204824366129e+3
    2.063468436776773e-001
    8.048303088078329e-002
   -3.667009475726021e-004
   -1.460201147413931e-003
   -2.586095375244759e-003
   -3.049813503085145e-007
    4.494611749252150e-005
    7.927512875033964e-005
   -1.235870224159925e-007
   -4.177551535814246e-009
   -4.302452311932423e-004
    6.337776244879493e-006
   -7.264046666691641e-010
   -5.107506824983828e-005
   -5.810472591789017e-009];

Idata = ~isnan(SP.*CT);
% ss = s(:); 
% ctt = ct(:);

A(:,1) = ones(length(SP(Idata)),1);
A(:,2) = CT(Idata);
A(:,3) = CT(Idata).^2;
A(:,4) = A(:,3).*CT(Idata);
A(:,5) = SP(Idata);
A(:,6) = SP(Idata).*CT(Idata);
A(:,7) = SP(Idata).*SP(Idata);

A(:,8) = CT(Idata);
A(:,9) = A(:,3);
A(:,10) = A(:,4);
A(:,11) = A(:,4).*CT(Idata);
A(:,12) = SP(Idata);
A(:,13) = A(:,6);
A(:,14) = SP(Idata).*A(:,4);
A(:,15) = SP(Idata).*sqrt(SP(Idata));
A(:,16) = A(:,15).*A(:,3);

gnum = A(:,1:7)*gcoeffs(1:7);

gden = 1+A(:,8:16)*gcoeffs(8:16);

g = nan(size(SP));
g(Idata) = gnum./gden; 

neutral_density_rf = g - 1000;


end
