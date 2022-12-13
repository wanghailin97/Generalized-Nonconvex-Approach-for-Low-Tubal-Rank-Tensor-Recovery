addpath(genpath(cd))
clear

pic_name = ['K:\WTNN\WTNN\LibADMM-master\image\testimg.jpg'];%输入你当前的文件名
I = double(imread(pic_name));
X = I/255;
    
[n1,n2,n3] = size(X);

opts.mu = 1e-3;
opts.tol = 1e-6;
opts.rho = 1.2;
opts.max_iter = 500;
opts.DEBUG = 1;


p = 0.3;
maxP = max(abs(X(:)));
omega = find(rand(n1*n2*n3,1)<p);
M = zeros(n1,n2,n3);
M(omega) = X(omega);


%% %% test TTNN,测试截断核范数
name='截断核范数';
disp(name);
r=3;
WTT=ones(n3);
alpha=0;
sign=-0.02;%对结果没有影响
[Xhat,obj,err,iter] = lrtc_wtnn(M,omega,r,WTT,alpha,sign,opts);

err
iter
obj
Xhat = max(Xhat,0);
Xhat = min(Xhat,maxP);
RSE = norm(X(:)-Xhat(:))/norm(X(:))
psnr = PSNR(X,Xhat,maxP)

figure(1)
subplot(1,3,1)
imshow(X/maxP)
subplot(1,3,2)
imshow(M/maxP)
subplot(1,3,3)
imshow(Xhat/maxP)


 pause

 %%测试加权核范数  WTNN
% % test lrtcR_snn
name='加权核范数';
disp(name);
r=0;
[U,S,V]=svd(M(:,:,1));
WTT=diag(S);%WTT可以输入其他不同值
alpha=1;
sign=-0.02;%对结果没有影响
[Xhat,obj,err,iter] = lrtc_wtnn(M,omega,r,WTT,alpha,sign,opts);

err
iter
obj

Xhat = max(Xhat,0);
Xhat = min(Xhat,maxP);
RSE = norm(X(:)-Xhat(:))/norm(X(:))
psnr = PSNR(X,Xhat,maxP)

figure(1)
subplot(1,3,1)
imshow(X/maxP)
subplot(1,3,2)
imshow(M/maxP)
subplot(1,3,3)
imshow(Xhat/maxP)
