function [MatrixL_3] = LinearL(j)
% Form the linear map from n^3 -> n^2 (on page 19)
n=2*j+1;
MatrixL_3=zeros(n^2,n^3);

for i1=1:1:n^3
    % we only concern about the first three dimension. So the index range are 1-n^3 corresponding to [1,1,1]-[n,n,n]
    [x,y,z]=ind2sub([n,n,n],i1);  
    
    % for any V[x,y,z,0,0,0], we can calculate the corresponding [m,l,m']
    m=j-x+1;
    l=-(j-y+1);
    m1=j-z+1;
    
    % Function "ComputeCo(j,m,m')" output the coefficient vector of (m,m')entry in the Wigner D-matrix
    coeffi=ComputeCo(j,m,m1); 
    
    % for given [m,l,m'], we can obtain the coefficient el
    el=coeffi(l+j+1);             
        
    % Find out the location in matrix L corresponding to ([m,m'],[x,y,z,0,0,0])
    colindex_3=sub2ind([n,n,n],x,y,z);    
    rowindex=sub2ind([n,n],m+j+1,m1+j+1);    
    MatrixL_3(rowindex,colindex_3)=el;
    
end