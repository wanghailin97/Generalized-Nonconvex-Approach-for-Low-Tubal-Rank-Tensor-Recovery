function [X,obj,err,iter] = lrtc_wtnn(M,omega,r,WTT,alpha,sign,opts)
%一般加权方法，初始权重使用第一次奇异值分解的倒数
% Solve the Low-Rank Tensor Completion (LRTC) based on Weighted Tensor Nuclear Norm (WTNN) problem by M-ADMM
%Input the weight parameter. This is the fixed weight parameter method.
%
%the weight parameter do not epsilon, the sign=0, else alpha=negative,sign=1
% min_X ||X||_{\omega}, s.t. P_Omega(X) = P_Omega(M)
%
% ---------------------------------------------
% Input:
%              M       -    d1*d2*d3 tensor
%            omega   -    index of the observed entries
%             WTT        -
%             是权重因子，是一个向量,可以根据实际情况手动输入，或者在代码中更换相关行，变成自动跟新权重的方法
%            apha        -   取负数，决定epsilon的大小
%            sign        -   决定是否使用epsilon
%              r          -   若其不为0，则表示矩阵的截断秩，即TTNNr，若其为0，按一般加权核范数计算
%             opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1

% Output:
%       X       -    d1*d2*d3 tensor
%       err     -    residual
%       obj     -    objective function value
%       iter    -    number of iterations

% version 1.0 - 18/09/2019
%
% Written by Canyi Lu (canyilu@gmail.com)
% The code is modificated by Mu(mylyw1314@tju.edu.cn)
% If you use this code, please cite the


% [n1,n2,n3]=size(M);
tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;

if ~exist('opts', 'var')
    opts = [];
end

if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end




dim = size(M);
k = length(dim);
omegac = setdiff(1:prod(dim),omega);
modelP=length(omega);
modelPM=modelP/(dim(1)*dim(2)*dim(3));%当使用比三阶更高的张量时，需要修改
epsilon=sign*exp(alpha*modelPM);
% epsilon=1e-8;
X = zeros(dim);
X(omega) = M(omega);

E = zeros(dim);
Y = E;


iter = 0;

for iter = 1 : max_iter
    Xk = X;
    Ek = E;
    % update X
    if r==0
    %% 不更新权重的情形
     [X,tnnX,~]=gprox_wtnn(-E+M+Y/mu,WTT,1/mu,epsilon); 
%% 更新权重的情形   
% [X,tnnX,WTT] = gprox_wtnn(-E+M+Y/mu,WTT,1/mu,epsilon);% 每次迭代加权
%        WT=WTT;
    else
        [X, tnnX] = gprox_ttnn(-E+M+Y/mu,1/mu,r);
    end
  %% update E
    E = M-X+Y/mu;
    E(omega) = 0;
 
    dY = M-X-E;    
    chgX= max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chg = max([chgX chgE max(abs(dY(:)))]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnX;
            err = norm(dY(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);    
end
obj = tnnX;
err = norm(dY(:));

 