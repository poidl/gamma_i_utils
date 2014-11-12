function [irow,jcol,n_bdy,bdy]=matrix_ij_bdy(wet,ibb)

[nz,ny,nx]=size(wet);
wet=wet(:);

sreg=cumsum(wet);
sreg(~wet)=nan;

% boundary
bdy=false(nz,ny,nx);
bdy(:,ibb)=true;
%keyboard
bdy= wet & bdy(:);

jcol= sreg(bdy); % column indices for matrix coef. 1
irow=(1:sum(bdy))';
n_bdy=sum(bdy);
%keyboard



     
end    


