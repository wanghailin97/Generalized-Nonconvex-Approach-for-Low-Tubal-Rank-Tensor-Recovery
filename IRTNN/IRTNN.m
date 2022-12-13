function X = IRTNN(fun, y, M, sizeX, gamma, err, x_initial, normfac, insweep, tol, decfac)

%% parameter settings
hfun_sg = str2func( [ fun '_sg'] ); 

if nargin < 5;   gamma = 1;                       end
if nargin < 6;   err = 1e-6;                      end
if nargin < 7;   x_initial = zeros(prod(sizeX),1);  end
if nargin < 8;  normfac = 1;                     end
if nargin < 9;  insweep = 200;                   end
if nargin < 10;  tol = 1e-5;                      end
if nargin < 11;  decfac = 0.9;                    end
mu = 1.1*normfac;
x = x_initial;
lambdaInit = decfac*max(abs(M(y,2))); lambda = lambdaInit;
f_current = norm(y-M(x,1)) + lambda*norm(x,1);

%% main
while lambda > lambdaInit*tol
    for ins = 1:insweep    
        f_previous = f_current;
        x = x + (1/mu)*M(y - M(x,1),2);
        x_tensor = reshape(x,sizeX) ;       
        [X,~,~] = prox_wtnn(x_tensor,hfun_sg,gamma,lambda,1/mu) ;        
        x = X(:) ;        
        f_current = norm(y-M(x,1)) + lambda*norm(x,1) ;
        if norm(f_current-f_previous)/norm(f_current + f_previous)<tol
            break;
        end
    end
    if norm(y-M(x,1))<err
        break ;
    end
    lambda = decfac*lambda;
end



