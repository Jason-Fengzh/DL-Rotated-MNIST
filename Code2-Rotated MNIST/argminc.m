function [c_opt] = argminc(Y,A,alpha,beta,gamma,j)
% Argmin over c
% Computes the analytical solution

Z = A*WignerMatrixGlobal(alpha,beta,gamma,j);
Znorm = norm(Z);
rZ = real(Z);
rY = real(Y);
iZ = imag(Z);
iY = imag(Y);

x = (trace(rZ'*rY)+trace(iZ'*iY))/(Znorm^2);
y = (-trace(iZ'*rY)+trace(rZ'*iY))/(Znorm^2);

% Optimal value
c_opt = x + 1j*y;

end