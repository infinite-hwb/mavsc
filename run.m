clc
clear,close all

%% Add path
addpath(genpath(cd))
addpath('.\tool');

%% Load Data 
load UCF
truth = gnd;
k = max(gnd);
K = KERNEL;
n = size(gnd,1);

%% MAVSC
fprintf('\nCentroid multiview M-AVSC\n');
mu=60;
lambda1=0.5;

% Coefficient matrix
[C] = centroid_MGSSC(n,K,mu,lambda1); 

% Obtain joint affinity matrix
A = BuildAdjacency(thrC(C,1));        

% Clustering result
grp = SpectralClustering(A, max(truth));
grps = bestMap(truth,grp);

acc = compacc(grps',gnd')
