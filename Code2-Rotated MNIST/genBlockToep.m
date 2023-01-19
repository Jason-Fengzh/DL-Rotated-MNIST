function [T] = genBlockToep(alpha,beta,gamma,j)
% Generate a Block Toeplitz tensor
n = 2*j + 1;
T = zeros(n,n,n);

for i1 = 1 : n
    for i2 = 1 : n
        for i3 = 1 : n
            c = alpha*(i1-j-1) + beta*(i2-j-1) + gamma*(i3-j-1);
            T(i1,i2,i3) = exp(1i*c);
        end
    end
end

end

