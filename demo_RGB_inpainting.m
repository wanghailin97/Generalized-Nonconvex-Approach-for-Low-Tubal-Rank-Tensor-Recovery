%% Demo of RGB image inpainting %%
%  test for IRTNN based RGB image inpainting and comparison with other SOTA methods
%==========================================================================
% This script compares several very realted methods
% listed as follows:
%   1. SNN
%   2. TNN
%   3. WTNN
%   4. PSTNN
%   5. IRNN
%   5. IRTNN (proposed method)
%
% -- three quality assessment (QA) indices -- PSNR, SSIM, FSIM, 
% -- are calculated for each method's denoising result, see Img_QA
%
% You can:
%       1. Type 'Demo' to to run various methods and see the pre-computed results.
%       2. Change test Img by simply modifying variable 'dataNum' in Demo.m 
%          (NOTE: make sure your Img meets the format requirements).
%       3. Change missing level by modifying variables  'missing_rate' in Demo.m
%       4. Select competing methods by turn on/off the enable-bits in Demo.m
%
% More detail can be found in 
% Hailin Wang et.al, Generalized Nonconvex Approach for Low-Tubal-Rank Tensor Recovery, TNNLS, 2022
%
% by Hailin Wang, 2022
%==========================================================================
% The inpainting results of the "1" RGB image when missing_rate = 30%
% ================== QA Results =====================
%    Method     PSNR     SSIM     FSIM      Time  
%  Observed   11.933    0.095    0.513    0.000   
%       SNN   31.376    0.906    0.942    9.039   
%       TNN   32.285    0.909    0.949    194.124   
%      WTNN   32.266    0.905    0.948    7.076   
%     PSTNN   31.366    0.890    0.942    25.088   
%      IRNN   28.680    0.835    0.916    34.060   
%     IRTNN   33.577    0.915    0.957    29.696  
% =========================================================================
% The inpainting results of the "7" RGB image when missing_rate = 70%
% ================== QA Results =====================
%    Method     PSNR     SSIM     FSIM      Time  
%  Observed   11.547    0.071    0.353    0.000   
%       SNN   41.235    0.986    0.992    8.565   
%       TNN   45.642    0.993    0.996    189.685   
%      WTNN   44.016    0.989    0.995    8.274   
%     PSTNN   45.862    0.993    0.997    23.315   
%      IRNN   36.864    0.963    0.982    17.048   
%     IRTNN   47.599    0.994    0.997    20.791  
%==========================================================================
%% Demo start
clc;
clear;close all;
rng('default');
% rng(1997)
addpath(genpath('lib'));
addpath(genpath('dataset'));
dataNum = '1';  % number index of RGB image in file dataset
%  Please make sure the RGB image is a cubic of size [height, width, 3] and
%  in range [0, 1], if it is in [0, 255], one can divided by 255
%  You can use other tensor data such as hyperspectral image, video, face data, network traffic data for test. 
%  Note some parameters might need to be reset for better performance

dataRoad = ['dataset/', dataNum, '.jpg'];
resultsRoad = ['results/RGB_Img_Inpainting/results_for_', dataNum];
if ~exist(resultsRoad); mkdir(resultsRoad); end
%% Set enable bits
Run_SNN       = 0;  % 1 or 0 <=> run or not
Run_TNN       = 0;  
Run_WTNN      = 0;  
Run_PSTNN     = 0;  
Run_IRNN      = 1;  
Run_IRTNN     = 1;  % our method

get_Recovered_Image = 1;  % if save the recovered image or not

%% Load Data 
methodName = {'Observed', 'SNN', 'TNN', 'WTNN', 'PSTNN', 'IRNN', 'IRTNN'};
Mnum = length(methodName);
data = double(imread(dataRoad))/255;  % load data
[height, width, band] = size(data);
sizeX = [height, width, band];

%% Observation
i = 1;
sampling_rate  = 0.3; % sampling_rate

disp(['=== the sampling rate is ', num2str(sampling_rate), ' ===']);
saveRoad = ['results/RGB_Img_Inpainting/results_for_', dataNum, '/', erase(num2str(sampling_rate),'.')];
if ~exist(saveRoad); mkdir(saveRoad); end
if exist([saveRoad, '/QA_Results.mat']); load([saveRoad, '/QA_Results.mat']); end
if exist([saveRoad, '/Results.mat']); load([saveRoad, '/Results.mat']); end
if get_Recovered_Image
    if ~exist([saveRoad,'/Recovered_Image']); mkdir([saveRoad,'/Recovered_Image']); end
    imwrite(data, [saveRoad, '/Recovered_Image/Original.jpg']);
end

m          = round(prod(sizeX)*sampling_rate);
sort_dim   = randperm(prod(sizeX));
Omega      = sort_dim(1:m); % sampling pixels' index
Obs        = zeros(sizeX);
Obs(Omega) = data(Omega); % observed Img

Results{i} = Obs;
[PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
enList = 1;

%% Run SNN
i = i+1;
if Run_SNN
    addpath(genpath(['competing methods/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the HaLRTC.m or example.m in SNN file
    Omega_SNN = zeros(sizeX);
    Omega_SNN(Omega) = 1;
    Omega_SNN = logical(Omega_SNN);
    alpha = [1, 1, 1e-3];
    alpha = alpha / sum(alpha);
    maxIter = 500;
    epsilon = 1e-5;
    beta = 1e-2;

    tic;
    Results{i} = HaLRTC(data, Omega_SNN, alpha, beta, maxIter, epsilon);
    Time(i) = toc;
%     % or we can equally use the lrtc_snn.m via ADMM sovler from the LibADMM-master
%     addpath(genpath('competing methods/LibADMM-master'));
%     opts = [];
%     alpha = [15,15,1]; % recommend by Lu et al., 
%     Results{i} = =lrtc_snn(Obs,Omega,alpha,opts);
%     Time(i) = toc;
%     rmpath(genpath('competing methods/LibADMM-master'));
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/RGB_Img_Inpainting/', methodName{i}]));
    enList = [enList, i];
end

%% Run TNN
i = i+1;
if Run_TNN
    addpath(genpath(['competing methods/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the tensor_cpl_admm.m or demo_tensor_completion.m in TNN file

    Omega_TNN     =        zeros(sizeX)                  ;
    known         =        sort_dim(1:m)                 ;
    Omega_TNN(known) =     1                             ;
    Omega_TNN     =        (Omega_TNN > 0)               ;
    normalize     =        max(Obs(:))                   ;
    inX           =        Obs/normalize                 ;
    gamma         =        1                             ;
    maxItr        =        1000                          ; % maximum iteration
    rho           =        0.01                          ;
    myNorm        =        'tSVD_1'                      ; % dont change for now
    A             =        diag(sparse(double(Omega_TNN(:)))); % sampling operator
    b             =        A * inX(:)                    ; % available data
    
    tic
    Results{i}    =    tensor_cpl_admm( A , b , rho , gamma , sizeX , maxItr , myNorm , 1 );
    Results{i}    =    Results{i} * normalize    ;
    Results{i}    =    reshape(Results{i},sizeX) ;
    Time(i) = toc;
    
    % or we can equally use the lrtc_tnn.m via ADMM sovler from the LibADMM-master by Lu et al.,
    % we recommend using the lrtc_tnn since it is much faster 
%     addpath(genpath('competing methods/LibADMM-master'));
%     opts = [];
%     Results{i} =lrtc_tnn(Obs,Omega,opts);
%     Time(i) = toc;
%     rmpath(genpath('competing methods/LibADMM-master'));

    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/RGB_Img_Inpainting/', methodName{i}]));
    enList = [enList, i];
end

%% Run WTNN
i = i+1;
if Run_WTNN
    addpath(genpath(['competing methods/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the lrtc_wtnn.m or testwtnn.m in WTNN file
    opts = [];
    opts.mu = 1e-3;
    opts.tol = 1e-6;
    opts.rho = 1.2;
    opts.max_iter = 500;
    opts.DEBUG = 1;
    r=0;
    [U,S,V]=svd(Obs(:,:,1));
    WTT=diag(S); 
    alpha=1;
    sign=-0.02; 
    
    tic
    Results{i} = lrtc_wtnn(Obs,Omega,r,WTT,alpha,sign,opts);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/RGB_Img_Inpainting/', methodName{i}]));
    enList = [enList, i];
end

%% Run PSTNN
i = i+1;
if Run_PSTNN
    addpath(genpath(['competing methods/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the LRTC_PSTNN.m or Demo_realworlddata_completion.m in PSTNN file
    opts = [];
    opts.mu = 1e-3;
    opts.tol = 1e-7;
    opts.rho = 1.1;
    opts.max_iter = 500;
    opts.DEBUG = 0;
    opts.max_mu = 1e10;
    tic
    [rankN,~] = prox_rankN(data,0.01);
    Results{i} =  LRTC_PSTNN(Obs,Omega,opts,rankN);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/RGB_Img_Inpainting/', methodName{i}]));
    enList = [enList, i];
end

%% Run IRNN
i = i+1;
if Run_IRNN
    addpath(genpath(['competing methods/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the IRNN.m or main.m in IRNN file
    % since the IRNN can only handle the matrix data, here, we implement it on each band

    % choose the used penalty function in IRNN and corresponding para
    fun = 'lp';         gamma = 0.5 ;
    % fun = 'scad';       gamma = 100;
    % fun = 'logarithm';  gamma = 10;
    % fun = 'mcp';        gamma = 10;
    % fun = 'cappedl1';   gamma = 1000;
    % fun = 'etp';        gamma = 0.1;
    % fun = 'geman';      gamma = 10;
    % fun = 'laplace';    gamma = 10;
    
    tic
    Results{i} = zeros(sizeX);
    for j = 1:band
        omega = find(rand(height*width,1)<sampling_rate);
        M = opRestriction(height*width, omega);
        x = data(:,:,j);
        y = M(x(:),1);
        Results{i}(:,:,j) = IRNN(fun,y,M,height,width,gamma) ;
    end
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/RGB_Img_Inpainting/', methodName{i}]));
    enList = [enList, i];
end

%% Run IRTNN
i = i+1;
if Run_IRTNN
    addpath(genpath('IRTNN'));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the IRTNN.m or main.m in IRTNN file
    % the IRTNN is a tensor based generalized nonconvex approach compared with IRNN, especially, proved with global convergence rigorously
    % the IRTNN is far better than the IRNN in process tensor data
    
    % choose the used penalty function in IRTNN and corresponding para
    fun = 'lp';         gamma = 0.5;
    % fun = 'scad';       gamma = 100;
%     fun = 'logarithm';  gamma = 10;
    % fun = 'mcp';        gamma = 100;
%     fun = 'etp';        gamma = 0.1;
    % fun = 'cappedl1';   gamma = 1000;
    % fun = 'geman';      gamma = 10;
    % fun = 'laplace';    gamma = 10;
    
    omega = find(rand(prod(sizeX),1)<sampling_rate);
    M = opRestriction(prod(sizeX), omega);
    y = M(data(:),1);
    
    tic
    Results{i} = IRTNN(fun,y,M,sizeX,gamma) ;
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/RGB_Img_Inpainting/', methodName{i}]));
    enList = [enList, i];
end

%% Show result
fprintf('\n');    

% enList = [1:7];

fprintf('================== QA Results =====================\n');
fprintf(' %8.8s    %5.5s    %5.5s    %5.5s     %5.5s  \n',...
    'Method', 'PSNR', 'SSIM', 'FSIM',  'Time');
% enList = 1:Mnum;
for i = 1:length(enList)
    fprintf(' %8.8s   %5.3f    %5.3f    %5.3f    %5.3f   \n',...
        methodName{enList(i)}, PSNR(enList(i)), SSIM(enList(i)), FSIM(enList(i)), Time(enList(i)));
end

fprintf('================== Show Results =====================\n');
figure;
plt_w = 3;
plt_h = ceil((length(enList)+1)/plt_w);
subplot(plt_h, plt_w,1); imshow(data);title('Original');
for ii = 1:length(enList)
    hold on;subplot(plt_h, plt_w,ii+1); imshow(Results{enList(ii)});title([methodName{enList(ii)},' psnr:',num2str(PSNR(enList(ii)),4)]);
end

%% save results
All = [PSNR; SSIM; FSIM; Time];
save([saveRoad,'/QA_Results'], 'All', 'PSNR', 'SSIM', 'FSIM', 'Time');
save([saveRoad,'/Results'], 'Results');