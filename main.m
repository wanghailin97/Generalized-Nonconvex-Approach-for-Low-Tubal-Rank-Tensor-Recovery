%% test on the tensor completion problem
clc;
clear;close all;
rng('default');rng(1997);
addpath(genpath('IRTNN'));

%% simulated data generate
newDataFlag = 1;
if newDataFlag == 1    
    n = 40; n1 = n; n2 = n; n3 = 5;
    r = 5; % tubal rank
    ML = randn(n1,r,n3)/n1; MR = randn(r,n2,n3)/n2;
    %dr = (n1+n2-r)*r*n3;
    %m = 3*dr;
    %p = m/(n1*n2*n3);
end

%% observation
p = 0.5; % sampling rate
omega = find(rand(n1*n2*n3,1)<p); 
IDX = omega ;
sizeX = [n1,n2,n3] ;
M = opRestriction(prod(sizeX), IDX);
X = tprod(ML,MR);
x = X(:); % vectorization of original tensor
y = M(x,1); % vectorization of observed tensor

%% select the used nonconvex penalty function and corresponding para 
% sometimes one need to adjust the 'gamma' value for better performance
fun = 'lp';         gamma = 0.5;
% fun = 'scad';       gamma = 100;
% fun = 'logarithm';  gamma = 10;
% fun = 'mcp';        gamma = 10;
% fun = 'etp';        gamma = 0.1;
% fun = 'cappedl1';   gamma = 1000;
% fun = 'geman';      gamma = 10;
% fun = 'laplace';    gamma = 10;

%% run IRTNN
tic
XRec = IRTNN(fun,y,M,sizeX,gamma) ;
time = toc

RelErr = norm(X(:)-XRec(:),'fro')/norm(X(:),'fro') % relative error, which is very small, meaning exact recovery

tubalrank(XRec)





