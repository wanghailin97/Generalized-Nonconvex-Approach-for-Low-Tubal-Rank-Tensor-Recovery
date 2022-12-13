%--------------Brief description-------------------------------------------
% This demo contains the implementation of the experiment for ¡®¡¯Sensitivity
% to initialization¡®¡¯
% More details in:
% Tai-Xiang Jiang, Ting-Zhu Huang, Xi-Le Zhao, Liang-Jian Deng;
% ''A novel nonconvex approach to recover the low-tubal-rank tensor data:
% when t-SVD meets PSSV'' submitted to Applied Mathematical Modelling
% (AMM).
% Contact: taixiangjiang@gmail.com
% Date: 7th Feb. 2018

% The authors would like to express their sincere thanks to 
% Dr. Canyi Lu (https://sites.google.com/site/canyilu/) 
% and 
% Zemin Zhang (https://sites.google.com/site/jamiezeminzhang/) 
% for their generous sharing of their code.Meanwhile, we promiose that 
% our code will all be updated as .m format as soon as the publication of our paper.


addpath(genpath(cd))
clear
n = 25;
n3 = 20;
times=1000;
alpha = [1 1 1];
alpha = alpha/sum(alpha);

r=5;
rhos = 0.1;
opts.sratio = r;
P = normrnd(0,1/sqrt(n),int8([n,r,n3]));
Q = normrnd(0,1/sqrt(n),int8([r,n,n3]));
L0 = tprod(P,Q);
[n1,n2,n3] = size(L0);


Ln = zeros(n1,n2,n3);
ind = find(rand(n1*n2*n3,1)<0.9);
Ln(ind) = L0(ind);

opts.mu = 10^-3;
opts.tol = 1e-14;
opts.rho = 1.2;
opts.max_iter = 200;
opts.DEBUG = 0;
opts.max_mu = 1e10;

tic
for flag=1:1000
    rng(flag)
    rankN = 5*ones(1,n3);
    [Lhat2,~,~,~] = LRTC_PSTNN(Ln,ind,opts,rankN);%LRTC_PSTNN(Ln,ind,rankN,mu,max_mu,tol,rho,max_iter,DEBUG,flag+1000);
    nmrse2(flag) = sqrt( norm(Lhat2(:)-L0(:)) )/ sqrt( norm(L0(:)) );
    
    if flag == 1 || mod(flag, 50) == 0
        disp(['iteration ' num2str(flag) 'completed']);
    end
end
hist(nmrse2,100)
toc