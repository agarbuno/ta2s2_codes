% Preparing the worskpace, adding possible directories I might need to
% compare the solutions of AIMS with the solutions of other non stochastic
% methods like Nelder Mead. 

function [X,y] = initialize_ex(n, func)
% clear;

addpath('~/Dropbox/Phd/Research/Matlab/keane-codes/')
addpath('~/Dropbox/Phd/Research/Matlab/keane-codes/Sampling Plans/')
addpath('~/Dropbox/Phd/Research/Matlab/keane-codes/Constructing a Surrogate/')
addpath('~/Dropbox/Phd/Research/Matlab/keane-codes/Exploring and Exploiting a Surrogate/')

% Set seed for pseudo random replication
% rng(201054345, 'twister');
% rng default
% rng(108727, 'twister');
% rng(09051991, 'twister');

% Return to original matlab format
format

%% To recover points from a latin hypercube

if strcmp(func, 'wingweight') 
    k = 10;
else
    k = 2;
end

% Number of sample points
%n = 5;
% figure(1); clf;
% load('data.mat')
% Create sampling plan
% X = bestlh(n,k,50,20);
X = lhsdesign(n,k);
%ModelInfo.X = bestlh(n,k,50,20);

% Calculate observed data
for i = 1:n
    y(i,1) = feval(func, (X(i,:)));
%ModelInfo.y(i,1) = branin(ModelInfo.X(i,:));
end


%% To recover specific data points (control and reference)

% load('data.mat')

%% To read from a benchmark dataset

% SARCOS DATASET
% load('sarcos_inv.m')
% X = sarcos_inv(:, 1:21);
% y = sarcos_inv(:, 22);
% 
% index = randsample(size(X,2), m);
% X = X(index,:);
% y = y(index,:);

% GEMSA DATASET
% infile = '~/Dropbox/Datafun/model2/training.txt';
% outfile = '~/Dropbox/Datafun/model2/out_5_100.txt';
% X = dlmread(infile);
% X = bsxfun(@rdivide, bsxfun(@minus, X, min(X)), (max(X) - min(X)));
% y = dlmread(outfile);

%% One dimensional toy examples

% % Create sampling plan
% X = linspace(-4.35,4.35,n)';
% % X = -5 + 10 * lhsdesign(n,1);
% % Calculate observed data
% y = oakley(X); 

%
% ... standarize a matrix
%

