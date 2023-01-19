function [LN] = MapL_beta(N)
% This is a simpler version of MapL(N)
% Here we only map V(beta) nx1 to d_mm1(beta) nxn actually n^2 x 1
n = 2*N+1;
LN_row_number = (N+1)*(2*N+1)*(2*N+3)/3;
LN = zeros(LN_row_number,n);

for j=0:1:N
    
    ni = 2*j+1;
    Lj = zeros(ni^2,n);
    diff = N-j;
    for i1=1:1:ni^2
        [m,m1]=ind2sub([ni,ni],i1); 
%         ComputeCo(j,m,m1).'
        Lj(i1,1+diff:n-diff) = ComputeCo(j,m-j-1,m1-j-1).';
    end
    start_row = j*(2*j-1)*(2*j+1)/3+1;
    end_row = (j+1)*(2*j+1)*(2*j+3)/3;
    LN(start_row:end_row,:) = Lj;
    
end

end