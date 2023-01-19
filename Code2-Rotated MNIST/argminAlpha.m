function [optimalvalue] = argminAlpha(Y,A,c,beta,gamma,j)
% This function performs a grid search over [0,2pi) to obtain the optimal
% alpha
% Method is a binary search with two tiers

% Set the number of grid points
n_gridpts = 10;

% First pass, search alpha over n_gridpts
step1 = 2*pi / n_gridpts;
optimalvalue = 0;

distance = norm(Y-c*A*WignerMatrixGlobal(0,beta,gamma,j));
for k = 1 : 1 : n_gridpts
    alpha = step1*k;
    distancek = norm(Y-c*A*WignerMatrixGlobal(alpha,beta,gamma,j));
    if distancek < distance
        optimalvalue = alpha;
        distance = distancek;
    end
end

% Second pass, search over 2*n_gridpts
step2 = step1 / n_gridpts;
NewSearchStart = optimalvalue - step1;
for k = 1 : 1 : 2*n_gridpts
    alpha = NewSearchStart + step2*k;
    distancek = norm(Y-c*A*WignerMatrixGlobal(alpha,beta,gamma,j));
    if distancek < distance
        optimalvalue = alpha;
        distance = distancek;
    end
end

% Compute the optimal value modulo [0,2*pi)
integ = floor(optimalvalue / (2*pi));
fract = (optimalvalue / (2*pi)) - integ;
optimalvalue = fract * (2*pi);

end