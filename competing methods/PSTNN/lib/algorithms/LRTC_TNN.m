function [X,obj,err,iter] = LRTC_TNN(M,omega,opts,M_true)

% Solve the Low-Rank Tensor Completion (LRTC) based on Tensor Nuclear Norm (TNN) problem by ADMM

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');             tol            = opts.tol;                 end
if isfield(opts, 'max_iter');    max_iter   = opts.max_iter;        end
if isfield(opts, 'rho');            rho          = opts.rho;                end
if isfield(opts, 'mu');            mu           = opts.mu;                end
if isfield(opts, 'max_mu');    max_mu   = opts.max_mu;        end
if isfield(opts, 'DEBUG');      DEBUG     = opts.DEBUG;          end

dim = size(M);


X = rand(dim);
X(omega) = M(omega);
Y = X;
Lambda = zeros(dim);

iter = 0;
for iter = 1 : max_iter
    Xk = X;
    Yk = Y;
    % update X
    [X,tnnX] = prox_tnn(Y-Lambda/mu,1/mu); 
    % update Y
    Y = 1/mu*(mu*X+Lambda);
    Y(omega) = M(omega);
 
    dLam = X-Y;    
    chgX = max(abs(Xk(:)-X(:)));
    chgY = max(abs(Yk(:)-Y(:)));
    chg = max([chgX chgY max(abs(dLam(:)))]);
    time = toc;
    if DEBUG
        if iter == 1 || mod(iter, 2) == 0
            obj = tnnX;
            err = norm(dLam(:));
           disp(['iter ' num2str(iter) ' , mu=' num2str(mu)  ' , psnr= ' num2str(psnr(X, M_true))...
                ', obj=' num2str(obj) ', err=' num2str(err) ',time=' num2str(time)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Lambda = Lambda + mu*dLam;
    mu = min(rho*mu,max_mu);    
end
obj = tnnX;
err = norm(dLam(:));

 