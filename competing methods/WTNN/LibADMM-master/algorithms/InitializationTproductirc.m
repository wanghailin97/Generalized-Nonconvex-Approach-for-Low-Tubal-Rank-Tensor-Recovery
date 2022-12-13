function  [C1,C,omega]=InitializationTproductirc(n1,n3,rou,lambda)
 r=rou*n1;
 %r=0.2*n1;
%r=0.3*n1;
% r=0.4*n1;
sigma=1/sqrt(n1);
M=zeros(n1,n1,n3);

X1=random('norm',0,sigma,n1,r,n3);
X2=random('norm',0,sigma,r,n1,n3);
 C=tproduct(X1,X2);
 %m=0.1*n1*n1*n3;
 % m=0.2*n1*n1*n3;
  m=lambda*n1*n1*n3;
%  m=0.4*n1*n1*n3;
 omg=rand(size(C));
 
 t=m/numel(C);
 omega = find(omg>t);
 M(omg>t)=C(omg>t);
 C1=M;