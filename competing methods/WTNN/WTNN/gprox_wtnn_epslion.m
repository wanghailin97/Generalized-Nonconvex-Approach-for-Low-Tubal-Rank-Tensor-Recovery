function [X, tnn,WTT ,trank] = gprox_wtnn_epslion(Y,WT,rho)

% The proximal operator of the tensor nuclear norm of a 3 way tensor
%Without consider the epsilon in the iteration, which will make a bad
%result for the recovery.

% min_X rho*||X||_*+0.5*||X-Y||_F^2
%
% Y     -    n1*n2*n3 tensor
%
% X     -    n1*n2*n3 tensor
% tnn   -    tensor nuclear norm of X
% trank -    tensor tubal rank of X
%wt     - the weighted of the tensor
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

[n1,n2,n3] = size(Y);
n12 = min(n1,n2);
cn3=ceil((n3+1)/2);
Y = fft(Y,[],3);
U = zeros(n1,n12,n3);
V = zeros(n2,n12,n3);
S = zeros(n12,n12,n3);
trank = 0;
for i = 1 : cn3
   
    [U(:,:,i),s,V(:,:,i)] = svd(Y(:,:,i),'econ');
 
%     wt=diag(wt);
%     WT(:,:,i)=s;
    s = diag(s);
   wt=WT(1:length(s));
%    s= max(s-rho./(WT+epsilon),1); 
  s= max(s-rho./wt,1);
%   disp(s);
    S(:,:,i) = diag(s);
    tranki = length(find(s~=0));
    trank = max(tranki,trank);
end
for i=cn3+1:n3
    U(:,:,i)=conj(U(:,:,n3-i+2));
    V(:,:,i)=conj(V(:,:,n3-i+2));
    S(:,:,i)=S(:,:,n3-i+2);
end
U = U(:,1:trank,:);                                                                                                                                                
V = V(:,1:trank,:);
S = S(1:trank,1:trank,:);

% my
% X=U*S*V';
% X=iff(X,[],3);
U = ifft(U,[],3);
S = ifft(S,[],3);
V = ifft(V,[],3);

X = tprod( tprod(U,S), tran(V));
WTT=WT;
S = S(:,:,1);
tnn = sum(S(:)); % return the tensor nuclear norm of X
