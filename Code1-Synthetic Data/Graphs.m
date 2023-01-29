%% Two methods together with comparision
clear;clc;  % Reset
rng(3)  % Seed
j=1;
n=2*j+1;
c1=1;  % c1,c2 is the number of dictionary elements (atoms)
c2=5;
K=50;  % Number of data
N=10;  % Divide the interval to N1/N2 step to search the minimum
times=3;  % Times we repeat the binary search type of procedure in function 'infidist'
lmd=0.1;  % Penalty factor
iteration=20;  % Number of iteration
interval=[0,2*pi;0,2*pi;0,2*pi];  % Interval for alpha,beta,gamma

[Phi_true,ymatrix,angle] = Data(j,K);
% y is a nxnxK array, angle is a 3xK matrix

%% Vanilla DL with c1=1 atoms
y=zeros(n^2,K);
for i=1:1:K
    y(:,i) = reshape(ymatrix(:,:,i),n^2,1);
end

% Generate Phi randomly for initialization
Phi=rand(n^2,c1)+1j*rand(n^2,c1);
Phi=normalize(Phi);

% Record the distance in every iteration
d0=zeros(1,c1);
matrixPhi0=zeros(n,n,c1);
for k=1:1:c1
    matrixPhi0(:,:,k)=reshape(Phi(:,k),n,n);
    [d0(k),~] = infidist(Phi_true,matrixPhi0(:,:,k),N,interval,times);
end
dmean0=mean(d0(:));
dmax0=max(d0(:));
dmin0=min(d0(:));

d=zeros(iteration,c1);
dmean1=zeros(1,iteration);
dmin1=zeros(1,iteration);
dmax1=zeros(1,iteration);
for it=1:1:iteration
    x=zeros(c1,K);
    for i=1:1:K
        x(:,i) = argminX(y(:,i),Phi,c1,lmd);
    end

    Phi=y/x;
    Phi=normalize(Phi);
    
    matrixPhi=zeros(n,n,c1);
    for k=1:1:c1
        matrixPhi(:,:,k)=reshape(Phi(:,k),n,n);
        [d(it,k),~] = infidist(Phi_true,matrixPhi(:,:,k),N,interval,times);    
    end
    dmean1(it)=mean(d(it,:));
    dmax1(it)=max(d(it,:));
    dmin1(it)=min(d(it,:));
end
d=[d0;d];
dmean1=[dmean0,dmean1];
dmin1=[dmin0,dmin1];
dmax1=[dmax0,dmax1];

%% Vanilla DL with c2=5 atoms
y=zeros(n^2,K);
for i=1:1:K
    y(:,i)=reshape(ymatrix(:,:,i),n^2,1);
end

% generate Phi randomly for initialization
Phi=rand(n^2,c2)+1j*rand(n^2,c2);
Phi=normalize(Phi);

d0=zeros(1,c2);
matrixPhi0=zeros(n,n,c2);
for k=1:1:c2
    matrixPhi0(:,:,k)=reshape(Phi(:,k),n,n);
    [d0(k),~] = infidist(Phi_true,matrixPhi0(:,:,k),N,interval,times);
end
dmean0=mean(d0(:));
dmax0=max(d0(:));
dmin0=min(d0(:));

d=zeros(iteration,c2);
dmean2=zeros(1,iteration);
dmin2=zeros(1,iteration);
dmax2=zeros(1,iteration);
for it=1:1:iteration
    x=zeros(c2,K);
    for i=1:1:K
        x(:,i) = argminX(y(:,i),Phi,c2,lmd);
    end

    Phi=y/x;
    Phi=normalize(Phi);
    
    matrixPhi=zeros(n,n,c2);
    for k=1:1:c2
        matrixPhi(:,:,k)=reshape(Phi(:,k),n,n);
        [d(it,k),~] = infidist(Phi_true,matrixPhi(:,:,k),N,interval,times);    
    end
    dmean2(it)=mean(d(it,:));
    dmax2(it)=max(d(it,:));
    dmin2(it)=min(d(it,:));
end
d=[d0;d];
dmean2=[dmean0,dmean2];
dmin2=[dmin0,dmin2];
dmax2=[dmax0,dmax2];

%% Dictionary learning with symmetries
X=randn(n);
Y=randn(n);
Z=X+1i*Y;
phi=Z/(norm(Z,'fro'));
% initialize a random normalized phi
z=zeros(n,n,K);
d1=zeros(1,iteration);
d2=ones(1,iteration);
d=zeros(1,iteration);
[d0,~]=infidist(Phi_true,phi,N,interval,times);
for m=1:1:iteration
    for i=1:1:K
        [z(:,:,i)] = argminiZ(ymatrix(:,:,i),phi,lmd);
    end
    % step1, update z
    [phi] = argminPhi(ymatrix,z);
    % step2, update phi
    [d1(m),~] = infidist(Phi_true,phi,N,interval,times);
    d(m)=min(d1(m),d2(m));
end
d=[d0,d];

%% draw the curves
x= 0:1:iteration;
errorbar(x,dmean1,dmean1-dmin1,dmax1-dmean1);
hold on
errorbar(x,dmean2,dmean2-dmin2,dmax2-dmean2,'--');
plot(x,d,'c-o');
xlabel('Iterate','FontSize',20);
ylabel('Distance','FontSize',20);
legend('Vanilla DL with 1 atom','Vanilla DL with 5 atoms','DL with symmetries (Our method)');
title('Irreducible of size 3\times3');
