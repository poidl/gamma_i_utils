load('data/lsqr_input.mat')


% A=[[1 1 1 0];...
%    [0 1 0 0];...
%    [0 0 1 0];...
%    [0 0 0 1]]; 
% A=sparse(A);

regions=get_coupled_regions(A);
keyboard