function [y] = Reshape_EmbedToAmbient(X,j)
% Maps (flattens) an input block diagonal matrix to the ambient space 

% Compute the ambient dimension
Fourier_Dim = (j+1)*(2*j+1)*(2*j+3)/3;

% Output
y = zeros(Fourier_Dim,1);

% Iteratively map
currambientloc = 0;
currembedloc = 0;

for ii = 0 : j
    % Get current block
    G = X(currembedloc+1:currembedloc+2*ii+1,currembedloc+1:currembedloc+2*ii+1);
    y(currambientloc+1:currambientloc+(2*ii+1)^2) = G(:);

    % Move dimensions
    currambientloc = currambientloc + (2*ii+1)^2;
    currembedloc = currembedloc + 2*ii + 1;
end

end