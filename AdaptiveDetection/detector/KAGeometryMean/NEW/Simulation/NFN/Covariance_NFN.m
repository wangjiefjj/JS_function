clc
clear 
close all
warning off
n = 0.5; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
tau_m = mu/(lambda-1);
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_tt = 0.9;
sigma_t = sqrt(sigma_tt);
rou = 0.90;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
Sigma = fun_rho(rou,N,1);
SNRout=-5:1:20; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% 系统导向矢量
R_KA = zeros(size(Sigma));
for i = 1:10000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + Sigma.*(t*t')/10000;
end
tic
RR = zeros(N,N);
parfor i =1:1e3
    Train = fun_TrainData(str_train,N,L,Sigma,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    [x0,tau] = fun_TrainData(str_train,N,1,Sigma,lambda,mu,opt_train); 
%     R_KA = zeros(size(Sigma));
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
%     R_KA = Sigma.*(t*t');
    R_KA2 = (Sigma).*(t*t');
%     R_KA = eye(N);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_LogM = fun_RLogEMean(Train,10);
    R_E = fun_RPowerEMean(Train,1,10);
    R_P = fun_RPowerEMean(Train,-1,10);
    [R_CC,alpha_cc(i)]=fun_CC(Train,fun_SCMN(Train),R_KA2);
    [R_LogCC,alpha_lecc(i)]=fun_LogCC_new(Train,R_KA,10);
    [R_ECC,alpha_ecc(i)]=fun_PowerCC(Train,R_KA,1,10);
    [R_PCC,alpha_pcc(i)]=fun_PowerCC(Train,R_KA,-1,10);
    R_SFP = fun_SFP(Train,1);
    error_E(i) = norm((R_E)-(Sigma),'fro');
    error_ECC(i) = norm((R_ECC)-(Sigma),'fro');
    error_P(i) = norm(R_P-Sigma,'fro');
    error_PCC(i) = norm(R_PCC-Sigma,'fro');
    error_LogCC(i) = norm(R_LogCC-Sigma,'fro');
    error_LogM(i) = norm((R_LogM)-(Sigma),'fro');
    
    error_RSFP(i) = norm(R_SFP-Sigma,'fro');
    error_RCC(i) = norm(R_CC-Sigma,'fro');
end
toc
m_errorE= mean(error_E)/norm(Sigma,'fro');
m_errorECC= mean(error_ECC)/norm(Sigma,'fro');
m_errorLogCC = mean(error_LogCC)/norm(Sigma,'fro');
m_errorP = mean(error_P)/norm(Sigma,'fro');
m_errorPCC = mean(error_PCC)/norm(Sigma,'fro');
m_errorLogM = mean(error_LogM)/norm(Sigma,'fro');
m_errorRSFP = mean(error_RSFP)/norm(Sigma,'fro');
m_errorRCC = mean(error_RCC)/norm(Sigma,'fro');

% 
% 
% mean_alpha_ecc = mean(alpha_ecc);
mean_alpha_lecc = mean(alpha_lecc);
mean_alpha_ecc = mean(alpha_ecc);
mean_alpha_pcc = mean(alpha_pcc);
mean_alpha_cc = mean(alpha_cc);

str = [str_train,'_Rerror_',num2str(n),'N','_s',num2str(sigma_tt),'.mat'];
save (str); 




