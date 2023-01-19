function [A_FD,A_Embed] = GenRandDict(j)
% Generate random dictionary
% A_Embed is a block diagonal matrix
% It acts by left multiplication on *matrices* with the same pattern
% A_FD is the SAME map, reshaped so that it acts by left multiplication
% on *vectors* where we remove the zero entries

% Generate the linear maps in embedded dimension and ambient dimension
A_Embed = zeros((j+1)^2,(j+1)^2);
Fourier_Dim = (j+1)*(2*j+1)*(2*j+3)/3;
A_FD = zeros(Fourier_Dim,Fourier_Dim);

c_dim = 0;
for ii = 0 : j
    % Generate each block randomly
    G = randn(2*ii+1,2*ii+1) + 1j*randn(2*ii+1,2*ii+1); % Random complex gaussian matrix
    A_Embed(ii^2+1:(ii+1)^2,ii^2+1:(ii+1)^2) = G;
    GK = kron(G,eye(2*ii+1)); % We take Kronecker products because we vectorize the argument
    % switched things around because we multiply on the right...
    A_FD(c_dim+1:c_dim+(2*ii+1)^2,c_dim+1:c_dim+(2*ii+1)^2) = GK;
    c_dim = c_dim + (2*ii+1)^2;
end

end

