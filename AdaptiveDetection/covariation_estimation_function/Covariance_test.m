clc
clear 
close all
warning off
n = 2; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 0.1;
rou = 0.9;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
Sigma = fun_rho(rou,N,1,0);
SNRout=-5:1:20; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
R_KA = zeros(size(Sigma));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + Sigma.*(t*t')/1000;
end
tic
opt = 4;
for i =1:1
    i
    Train = fun_TrainData(str_train,N,L,Sigma,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    [x0,tau] = fun_TrainData(str_train,N,1,Sigma,lambda,mu,opt_train); 
    %%协方差估计
    R_SCM = fun_SCMN(Train);
    
    R_LogM = fun_RLogEMean(Train,opt);
    R_LogCC = fun_LogCC_new(Train,R_KA,opt);
%     R_RPCC2 = fun_PowerCC(Train,R_KA,2,opt);
    R_RP1 = fun_RPowerEMean(Train,1,opt);%-1,0.5,2
    R_RP2 = fun_RPowerEMean(Train,2,opt);%-1,0.5,2
    R_FP = fun_FP(Train);
    [R_SFP1,rho] = fun_SFP(Train,1);
    R_SFP2 = fun_SFP(Train,2);
    R_NSCM = fun_NSCMN(Train);
    R_CC = fun_CC(Train,R_NSCM,R_KA);
    R_x0 = (fun_SCMN(x0));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_RSCM(i) = norm(R_SCM-Sigma,'fro');
    error_RLogM(i) = norm(R_LogM-Sigma,'fro');
    error_RLogCC(i) = norm(R_LogCC-Sigma,'fro');
    error_RP2(i) = norm(R_RP2-Sigma,'fro');
    error_RP1(i) = norm(R_RP1-Sigma,'fro');
%     error_RPCC2(i) = norm(R_RPCC2-Sigma,'fro');
    error_RFP(i) = norm(R_FP-Sigma,'fro');
    error_RSFP1(i) = norm(R_SFP1-Sigma,'fro');
    error_RSFP2(i) = norm(R_SFP2-Sigma,'fro');
    error_RNSCM(i) = norm(R_NSCM-Sigma,'fro');
end
toc
m_errorR_KA = norm(R_KA-Sigma,'fro')/norm(Sigma,'fro');
m_errorRSCM = mean(error_RSCM)/norm(Sigma,'fro');
m_errorRLogM = mean(error_RLogM)/norm(Sigma,'fro');
m_errorRLogCC = mean(error_RLogCC)/norm(Sigma,'fro');
m_errorRP2 = mean(error_RP2)/norm(Sigma,'fro');
m_errorRP1 = mean(error_RP1)/norm(Sigma,'fro');
% m_errorRPCC2 = mean(error_RPCC2)/norm(Sigma,'fro');
m_errorRFP = mean(error_RFP)/norm(Sigma,'fro');
m_errorRSFP1 = mean(error_RSFP1)/norm(Sigma,'fro');
m_errorRSFP2 = mean(error_RSFP2)/norm(Sigma,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(Sigma,'fro');


% num = 1:1000;
% plot(num,alpha_cc,'b')
% hold on
% plot(num,alpha_ML,'g')
% plot(num,alpha_ecc,'r')
% plot(num,alpha_lecc,'k')
% legend('CC','ML','E','LE')




