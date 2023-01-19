%% Script for SO3-invariant Dictionary Learning
% Experiments to learn some dictionary for digits

%% Load data

% Reset, in case we mix up results from previous runs
clear all; close all; clc;

% Change the data_name, 'digi' for i=0,1,...,9
load('fourier_coefficients_16.mat','dig9')

% This is the bandwidth
% Bandwidth of 6 + 1 coefficients
N = 10; % Can be 6 or less

% get data with size corresponding to N
Data = dig9(:,1:(N+1)^2,1:(N+1)^2);
% Set it to a common name that we use throughout the script. We only have to modify the filename here..

%% Re-shaping data, and other pre-processing steps

n_datapts = 200; % Using all data-points for now

nInnerLoopIterates = 5;

% Setting up the Coding step -- global scope

global ToepToWignerL_alpha
ToepToWignerL_alpha = MapL_alpha(N);

global ToepToWignerL_beta
ToepToWignerL_beta = MapL_beta(N);

global ToepToWignerL_gamma
ToepToWignerL_gamma = MapL_gamma(N);

Fourier_Dim = (N+1)*(2*N+1)*(2*N+3)/3;
Y = zeros(Fourier_Dim,n_datapts);
Yt = zeros(Fourier_Dim,n_datapts);
for ii = 1 : n_datapts
    f = reshape(Data(ii,:,:),(N+1)^2,(N+1)^2);
    Y(:,ii) = Reshape_EmbedToAmbient(f,N);
    Yt(:,ii) = Reshape_EmbedToAmbient(f',N); % Creating a transposed version
end

%% Initialize Dictionary

[A_FD,A_Embed] = GenRandDict(N);

%% The Dictionary Learning Procedure

n_Iterates = 10;
n_progressbar = floor(n_datapts/20); % The progess bar has 20 dots

fprintf('|It#|....................|Sq.Err|Post.Norm.Err|\n')

err_ii = 0; % break condition

for ii = 1 : n_Iterates
    tic
    fprintf('|#');
    fprintf(int2str(ii))
    fprintf('|')

    % The Coding Step
    X_Ambient = zeros(Fourier_Dim,n_datapts);
    X_Embed = zeros(n_datapts,(N+1)^2,(N+1)^2);
    for jj = 1 : n_datapts
        X_Ambient(:,jj) = Coding(Yt(:,jj),A_FD,N,nInnerLoopIterates); % Switch out to Yt
        X_Embed(jj,:,:) = Reshape_AmbientToEmbed(X_Ambient(:,jj),N);
        if mod(jj,n_progressbar) == 0
            fprintf('.'); % This is the progress bar
        end
    end
    fprintf('|');
    
    err1 = norm(Yt - A_FD*X_Ambient,'fro')^2/norm(Y,'fro')^2;
    fprintf(num2str(err1))
    fprintf('|')
    
    % The Dictionary Update Step
    [A_FD,A_Embed] = DictUpdate(Yt,X_Embed,N);
    
    err3 = norm(Yt - A_FD*X_Ambient,'fro')^2/norm(Y,'fro')^2;
    fprintf(num2str(err3))
    fprintf('|')
    
    % Save output
    fprintf(' Dict update... ');
    fname = strcat("SO3-Dict_Ambient_it",int2str(ii),".mat");
    save(fname,'A_FD')
    fname = strcat("SO3-Dict_Embed_it",int2str(ii),".mat");
    save(fname,'A_Embed')

    fprintf(' Iter completed and saved.\n');
    
    toc
    
%     % set a break conditon
%     if abs(err_ii - err1) <= 0.000001
%         break;
%     end
%     err_ii = err1;
    
end

% % Save the data of the last Iteration
% fname = strcat("SO3-Dict_Ambient_it",int2str(ii),".mat");
% save(fname,'A_FD')
% fname = strcat("SO3-Dict_Embed_it",int2str(ii),".mat");
% save(fname,'A_Embed')
