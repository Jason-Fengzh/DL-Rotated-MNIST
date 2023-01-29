%% Script for Dictionary Learning with symmetries (our method)
clear;clc;  % Reset
% Initializtion
rng(3)  % Seed
j=2;
n=2*j+1;
K=50; % Number of data
N1=10;  % Divide the interval to N1/N2 step to search the minimum
N2=15; 
times=3; % For calculating the distance, times we repeat the binary search
lmd=0.1;  % Penalty factor
iteration=20;  % Number of iteration

%% Generate data
[phi_true,y,angle]= Data(j,K);
% y is nxnxK array, angle is 3xK matrix

% initialize a random normalized phi
X=randn(n);
Y=randn(n);
Z=X+1i*Y;
phi=Z/(norm(Z,'fro'));

%% Compute the distance betweem true value and experiment result
z=zeros(n,n,K);
d1=zeros(1,iteration);
d2=ones(1,iteration);
d=zeros(1,iteration);
interval=[0,2*pi;0,2*pi;0,2*pi];
[d0,~]=infidist(phi_true,phi,N1,interval,times);
for m=1:1:iteration
    for i=1:1:K
        [z(:,:,i)] = argminiZ(y(:,:,i),phi,lmd);
    end
    % step1, update z
    [phi] = argminPhi(y,z);
    % step2, update phi
    [d1(m),~] = infidist(phi_true,phi,N1,interval,times);
%     if(d1(i)>0.05)
%         [d2(i),~] = infdist(y(:,:,i),phi,N2,interval,3);        
%     end
    d(m)=min(d1(m),d2(m));
end
d=[d0,d];

%% draw the curve and save the data
x= 0:1:iteration;
plot(x,d,'-o');
xlabel('Iterate','FontSize',20);
ylabel('Distance','FontSize',20);

% save('d_phi_y_z.mat','d','phi','phi_true','y','z');
