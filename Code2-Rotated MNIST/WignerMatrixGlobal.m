function [D] = WignerMatrixGlobal(alpha,beta,gamma,j)
% Generate the Wigner D matrix
% This current version outputs D as a tall vector with dimension (j+1)*(2*j+1)*(2*j+3)/3

n = 2*j + 1;
global ToepToWignerL_alpha
global ToepToWignerL_beta
global ToepToWignerL_gamma

Valpha = zeros(n,1);
Vbeta = zeros(n,1);
Vgamma = zeros(n,1);
for m=1:1:n
    Valpha(m)=exp(1i*alpha*(m-j-1));
    Vbeta(m)=exp(1i*beta*(m-j-1));
    Vgamma(m)=exp(1i*gamma*(m-j-1));
end


d2_vec = ToepToWignerL_beta * Vbeta;
multiplier_alpha = ToepToWignerL_alpha * Valpha;
multiplier_gamma = ToepToWignerL_gamma * Vgamma;

D = multiplier_alpha .* d2_vec .* multiplier_gamma;

end