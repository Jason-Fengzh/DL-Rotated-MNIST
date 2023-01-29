%% Script for Vanilla Dictionary Learning
clear;clc;  % Reset
% Initializtion
rng(3)  % Seed
j=2;
n=2*j+1;
c=5;  % Number of dictionary elements
K=50;  % Number of data
lmd=0.1;  % Penalty factor
iteration=20;  % Number of iteration
interval=[0,2*pi;0,2*pi;0,2*pi];  % Interval for alpha,beta,gamma
N=10;
times=3;

%% Generate data
[Phi_true,ymatrix,angle]= Data(j,K);

y=zeros(n^2,K);
for i=1:1:K
    y(:,i)=reshape(ymatrix(:,:,i),n^2,1);
end

% initialize a random normalized phi
Phi=rand(n^2,c)+1j*rand(n^2,c);
Phi=normalize(Phi);

%% Compute the distance betweem true value and experiment result
d0=zeros(1,c);
matrixPhi0=zeros(n,n,c);
for k=1:1:c
    matrixPhi0(:,:,k)=reshape(Phi(:,k),n,n);
    [d0(k),~] = infidist(Phi_true,matrixPhi0(:,:,k),N,interval,times);
end
dmean0=mean(d0(:));
dmax0=max(d0(:));
dmin0=min(d0(:));
% Compute and record in every iteration
d=zeros(iteration,c);
dmean=zeros(1,iteration);
dmin=zeros(1,iteration);
dmax=zeros(1,iteration);
for it=1:1:iteration
    x=zeros(c,K);
    for i=1:1:K
        x(:,i) = argminX(y(:,i),Phi,c,lmd);
    end

    Phi=y/x;
    Phi=normalize(Phi);    
    matrixPhi=zeros(n,n,c);
    for k=1:1:c
        matrixPhi(:,:,k)=reshape(Phi(:,k),n,n);
        [d(it,k),~] = infidist(Phi_true,matrixPhi(:,:,k),N,interval,times);    
    end
    dmean(it)=mean(d(it,:));
    dmax(it)=max(d(it,:));
    dmin(it)=min(d(it,:));
end
d=[d0;d];
dmean=[dmean0,dmean];
dmin=[dmin0,dmin];
dmax=[dmax0,dmax];

%% draw the curve and save the data
x= 0:1:iteration;
plot(x,dmean,'--');
errorbar(x,dmean,dmean-dmin,dmax-dmean,'--');
xlabel('Iterate','FontSize',20);
ylabel('Distance','FontSize',20);  

% save('d_phi_y.mat','dmean','dmax','dmin','matrixPhi','Phi_true','y');