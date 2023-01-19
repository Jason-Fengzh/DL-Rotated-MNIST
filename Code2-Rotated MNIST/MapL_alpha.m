function [LN] = MapL_alpha(N)
% A linear map from Valpha nx1 to the alpha_multiplier of a WignerMatrix
n = 2*N+1;
LN_row_number = (N+1)*(2*N+1)*(2*N+3)/3;
LN = zeros(LN_row_number,n);

for j=0:1:N
    
    ni = 2*j+1;
    Lj = zeros(ni^2,n);
    for i1=1:1:ni^2
        [m,~]=ind2sub([ni,ni],i1); 
        alpha_loc = N+1 - (m-j-1);
        Lj(i1,alpha_loc) = 1;
    end
    start_row = j*(2*j-1)*(2*j+1)/3+1;
    end_row = (j+1)*(2*j+1)*(2*j+3)/3;
    LN(start_row:end_row,:) = Lj;
    
end

end