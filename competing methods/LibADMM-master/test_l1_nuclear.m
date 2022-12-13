%
% References:
%
% C. Lu. A Library of ADMM for Sparse and Low-rank Optimization. National University of Singapore, June 2016.
% https://github.com/canyilu/LibADMM.
% C. Lu, J. Feng, S. Yan, Z. Lin. A Unified Alternating Direction Method of Multipliers by Majorization 
% Minimization. IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 40, pp. 527-541, 2018
%


addpath(genpath(cd))
clear

%% Examples for testing the sparse models
% For detailed description of the sparse models, please refer to the Manual.


%% generate toy data
d = 10;
na = 200;
nb = 100;

A = eye(d,na);

k = 0.2; % k-ratio spasity
omega = find(rand(na,nb)<k);
M = randn(na,nb);
X = zeros(na,nb);
X(omega) = M(omega);

B = A*X;
b = B(:,1);

opts.tol = 1e-6; 
opts.max_iter = 1000;
opts.rho = 1.1;
opts.mu = 1e-4;
opts.max_mu = 1e10;
opts.DEBUG = 0;

%% l1
[X2,obj,err,iter] = l1(A,B,opts);
figure,
stem(X2(:,1))
norm(X2(:)-X(:))

[X21,obj1,err1,iter1] = nuclear(A,B,opts);
figure,
stem(X21(:,1))
norm(X21(:)-X(:))

