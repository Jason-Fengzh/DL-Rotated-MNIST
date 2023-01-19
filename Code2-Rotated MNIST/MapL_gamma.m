function [LN] = MapL_gamma(N)
% A linear map from Vgamma nx1 to the gamma_multiplier of a WignerMatrix
n = 2*N+1;
LN_row_number = (N+1)*(2*N+1)*(2*N+3)/3;
LN = zeros(LN_row_number,n);

for j=0:1:N
    
    ni = 2*j+1;
    Lj = zeros(ni^2,n);
    for i1=1:1:ni^2
        [~,m1]=ind2sub([ni,ni],i1); 
        gamma_loc = N+1 - (m1-j-1);
        Lj(i1,gamma_loc) = 1;
    end
    start_row = j*(2*j-1)*(2*j+1)/3+1;
    end_row = (j+1)*(2*j+1)*(2*j+3)/3;
    LN(start_row:end_row,:) = Lj;
    
end

end