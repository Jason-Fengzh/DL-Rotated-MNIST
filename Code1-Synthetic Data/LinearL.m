%% The linear map from tensor space to Wigner D-matrix space
function [MatrixL_3] = LinearL(j)
n=2*j+1;
% For a map from n^3 -> n^2:
MatrixL_3=sparse(n^2,n^3);
for i1=1:1:n^3
    [x,y,z]=ind2sub([n,n,n],i1);  
    % we only concern about the first three dimension. So the index range are 1-n^3 corresponding to [1,1,1]-[n,n,n]
    m=j-x+1;
    l=-(j-y+1);
    m1=j-z+1;                    
    % for any V[x,y,z,0,0,0], we can calculate the corresponding [m,l,m']
    coeffi=ComputeCo(j,m,m1); 
    % function "ComputeCo(j,m,m')" output the coefficient vector of (m,m')entry in the Wigner D-matrix
    el=coeffi(l+j+1);             
    % for given [m,l,m'], we can obtain the coefficient el
    
    % Find out the location in matrix L corresponding to ([m,m'],[x,y,z,0,0,0])
    colindex_3=sub2ind([n,n,n],x,y,z);    
    rowindex=sub2ind([n,n],m+j+1,m1+j+1);       
    MatrixL_3(rowindex,colindex_3)=el;
    % MatrixL_3 is the map from n^3 -> n^2 we want
end