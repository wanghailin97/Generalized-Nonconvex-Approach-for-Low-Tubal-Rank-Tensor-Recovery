addpath(genpath(cd))
clear

pic_name = ['K:\WTNN\WTNN\LibADMM-master\image\testimg.jpg'];%�����㵱ǰ���ļ���
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


%% %% test TTNN,���ԽضϺ˷���
name='�ضϺ˷���';
disp(name);
r=3;
WTT=ones(n3);
alpha=0;
sign=-0.02;%�Խ��û��Ӱ��
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

 %%���Լ�Ȩ�˷���  WTNN
% % test lrtcR_snn
name='��Ȩ�˷���';
disp(name);
r=0;
[U,S,V]=svd(M(:,:,1));
WTT=diag(S);%WTT��������������ֵͬ
alpha=1;
sign=-0.02;%�Խ��û��Ӱ��
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
