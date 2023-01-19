function [A_FD,A_Embed] = DictUpdate(Y,X,j)
% Dictionary Update via Least Squares

% Generate the linear maps in embedded dimension and ambient dimension
Fourier_Dim = (j+1)*(2*j+1)*(2*j+3)/3;
A_FD = zeros(Fourier_Dim,Fourier_Dim);

[n,~,~] = size(X);

% Creating the Linear system
XX = zeros((j+1)^2,(j+1)^2);
YX = zeros((j+1)^2,(j+1)^2);

for ii = 1 : n
    x = reshape(X(ii,:,:),(j+1)^2,(j+1)^2);
    y = Reshape_AmbientToEmbed(Y(:,ii),j);
    XX = XX + x'*x;
    YX = YX + x'*y; % No transpose because the input data is Yt
end

% Solving the Linear System
A_Embed = (YX/XX)';

% Re-shape the solution, with Kronecker product, to the ambient space
c_dim = 0;
for ii = 0 : j
    G = A_Embed(ii^2+1:(ii+1)^2,ii^2+1:(ii+1)^2);
    GK = kron(G,eye(2*ii+1));
    A_FD(c_dim+1:c_dim+(2*ii+1)^2,c_dim+1:c_dim+(2*ii+1)^2) = GK;
    c_dim = c_dim + (2*ii+1)^2;
end

end

