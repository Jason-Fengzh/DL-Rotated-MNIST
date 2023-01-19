function [X] = Reshape_AmbientToEmbed(y,j)
% Maps (reshapes) an input vector in ambient space to the block diagonal
% matrix structure

% Output
X = zeros((j+1)^2,(j+1)^2);

% Iteratively map
currambientloc = 0;
currembedloc = 0;
for ii = 0 : j
    % Get current block
    g = y(currambientloc+1:currambientloc+(2*ii+1)^2);
    G = reshape(g,(2*ii+1),(2*ii+1));
    X(currembedloc+1:currembedloc+2*ii+1,currembedloc+1:currembedloc+2*ii+1) = G;

    % Move dimensions
    currambientloc = currambientloc + (2*ii+1)^2;
    currembedloc = currembedloc + 2*ii + 1;
end

end