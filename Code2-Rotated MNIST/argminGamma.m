function [optimalvalue] = argminGamma(Y,A,c,alpha,beta,j)
% This function performs a grid search over [0,2pi) to obtain the optimal
% gamma
% Method is a binary search with two tiers

% Set the number of grid points
n_gridpts = 10;

% First pass, search alpha over n_gridpts
step1 = 2*pi / n_gridpts;
optimalvalue = 0;

distance = norm(Y-c*A*WignerMatrixGlobal(alpha,beta,0,j));
for k = 1 : 1 : n_gridpts
    gamma = step1*k;
    distancek = norm(Y-c*A*WignerMatrixGlobal(alpha,beta,gamma,j));
    if distancek < distance
        optimalvalue = gamma;
        distance = distancek;
    end
end

% Second pass, search over 2*n_gridpts
step2 = step1 / n_gridpts;
NewSearchStart = optimalvalue - step1;
for k = 1 : 1 : 2*n_gridpts
    gamma = NewSearchStart + step2*k;
    distancek = norm(Y-c*A*WignerMatrixGlobal(alpha,beta,gamma,j));
    if distancek < distance
        optimalvalue = gamma;
        distance = distancek;
    end
end

% Compute the optimal value modulo [0,2*pi)
integ = floor(optimalvalue / (2*pi));
fract = (optimalvalue / (2*pi)) - integ;
optimalvalue = fract * (2*pi);

end