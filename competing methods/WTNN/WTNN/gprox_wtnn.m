function [X, tnn,WTT ,trank] = gprox_wtnn(Y,WT,rho,epsilon)

% The proximal operator of the tensor nuclear norm of a 3 way tensor
%This operator need input the weight parameter, give the next weight
%parameter for the next iteration. Consider the epsilon in the iteration.
%
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

[U(:,:,1),s,V(:,:,1)] = svd(Y(:,:,1),'econ');
WTT=diag(s);
wt=WT;
s = diag(s);
s= max(s-rho./(wt+epsilon),0);
S(:,:,1) = diag(s);
tranki = length(find(s~=0));
trank = max(tranki,trank);
for i = 2 : cn3  
    [U(:,:,i),s,V(:,:,i)] = svd(Y(:,:,i),'econ');
    s = diag(s);
%    s= max(s-rho./(WT+epsilon),1); 
 s= max(s-rho./(wt+epsilon),0);
 S(:,:,i) = diag(s);
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

S = S(:,:,1);
tnn = sum(S(:)); % return the tensor nuclear norm of X
