function [L,S,err,iter] = TRPCA_TNN(X,lambda,opts)

% Solve the Tensor Robust Principal Component Analysis based on Tensor Nuclear Norm problem by ADMM

% Many thanks for the Canyi Lu sharing the code, we update it at 7th Feb. 2018



if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
dim = size(X);
L = zeros(dim);
S = L;
Y = L;

iter = 0;
for iter = 1 : max_iter
    Lk = L;
    Sk = S;
    % update L
    [L,tnnL] = prox_tnn(-S+X-Y/mu,1/mu);
    % update S
    S = prox_l1(-L+X-Y/mu,lambda/mu);
  
    dY = L+S-X;
    chgL = max(abs(Lk(:)-L(:)));
    chgS = max(abs(Sk(:)-S(:)));
    chg = max([ chgL chgS max(abs(dY(:))) ]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnL+lambda*norm(S(:),1);
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
obj = tnnL+lambda*norm(S(:),1);
err = norm(dY(:));