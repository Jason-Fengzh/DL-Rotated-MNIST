function [X] = Coding(Y,A,j,nIterates)
% Sparse Coding Procedure
% We find the optimal c,alpha,beta,gamma minimizing
% ||Y - A c D(alpha,beta,gamma)||

% Number of iterates -- got to test more to find out what is a suitable
% choice here
%nIterates = 5;

% Initialize
alpha = 0.0;
beta = 0.0;
gamma = 0.0;
c = 0.0;

% Main coordinate descent procdure
for ii = 1 : nIterates
    alpha = argminAlpha(Y,A,c,beta,gamma,j);
    beta = argminBeta(Y,A,c,alpha,gamma,j);
    gamma = argminGamma(Y,A,c,alpha,beta,j);
    c = argminc(Y,A,alpha,beta,gamma,j);
end

% Final output
X = c*WignerMatrixGlobal(alpha,beta,gamma,j);

end