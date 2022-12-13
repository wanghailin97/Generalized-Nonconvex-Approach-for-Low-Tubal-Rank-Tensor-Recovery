%--------------Brief description-------------------------------------------
% This demo contains the implementation of the experiment for image recover
% More details in: 
% Tai-Xiang Jiang, Ting-Zhu Huang, Xi-Le Zhao, Liang-Jian Deng;
% ''A novel nonconvex approach to recover the low-tuba-rank tensor data:
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
clear all;clc;close all;

%% load the testing images
im_num = 1;%     1->starfish;	2->door;     3->hat1;	 4-> hat2

switch im_num
      case 1
            im_name = 'SwitchLight1049.bmp';
        case 2
            im_name = 'door.png';
        case 3
            im_name = 'hat1.png';
        case 4
            im_name = 'hat2.png';
end

X = im2double(imread(im_name));
if im_num>2   %cropping
    X = X(1:256,1:256,:);
end

[n1,n2,n3] = size(X);
maxP = max(abs(X(:)));
[n1,n2,n3] = size(X);
Xn = X;

rhos = 0.2;  %% The ratio of corrupted pixels

ind = find(rand(n1*n2*n3,1)<rhos);   
Xn(ind) = rand(length(ind),1);   %% corrupted by uniform distributed values

%% index of the corrupted image
PSNR0 = psnr(Xn,X,maxP);    
SSIM0 = ssim(Xn,X);

%% parameter setting
opts.mu = 1e-3;
opts.tol = 1e-7;
opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 0;
opts.max_mu = 1e10;

%% Tensor RRPCA based on SNN

alpha = [15 15 1.5];
[Xhat1,E1,~,~] = TRPCA_SNN(Xn,alpha,opts);

Xhat1 = max(Xhat1,0);
Xhat1 = min(Xhat1,maxP);
PSNR1 = psnr(Xhat1,X,maxP);
SSIM1 = ssim(Xhat1,X);


%% Tensor RRPCA based on TNN

[n1,n2,n3] = size(Xn);
lambda = 1/sqrt(max(n1,n2)*n3);
[Xhat2,E2,~,~] = TRPCA_TNN(Xn,lambda,opts);

 
Xhat2 = max(Xhat2,0);
Xhat2 = min(Xhat2,maxP);
PSNR2 = psnr(Xhat2,X,maxP);
SSIM2 = ssim(Xhat2,X);

%% Tensor RRPCA based on PSTNN

opts.rankN = prox_rankN1(X,0.02);
[n1,n2,n3] = size(Xn);
lambda = 1/sqrt(max(n1,n2)*n3);
[Xhat3,E3,~,~] = TRPCA_PSTNN(Xn,lambda,opts);

Xhat3 = max(Xhat3,0);
Xhat3 = min(Xhat3,maxP);
PSNR3 = psnr(Xhat3,X,maxP);
SSIM3 = ssim(Xhat3,X);

%% illustration of the results;
figure
subplot(2,3,1)
imshow(X);title('Original image');
subplot(2,3,2)
imshow(Xn);title('Observed image');
subplot(2,3,3)
imshow(Xhat3);title('Recovered by PSTNN');
% subplot(2,3,4)
% imshow(E3);title('The sparse corruption');
subplot(2,3,5)
imshow(Xhat1);title('Recovered by SNN');
subplot(2,3,6)
imshow(Xhat2);title('Recovered by TNN');

disp(['image name : ' im_name]);
disp(['Index || observed ||   SNN  ||   TNN  ||   PSTNN ']);
disp(['PSNR || ' num2str(PSNR0) ' ||  ' num2str(PSNR1) ' ||  ' num2str(PSNR2) ' ||  ' num2str(PSNR3) ]);
disp(['SSIM || ' num2str(SSIM0) ' ||  ' num2str(SSIM1) ' ||  ' num2str(SSIM2) ' ||  ' num2str(SSIM3)  ]);















