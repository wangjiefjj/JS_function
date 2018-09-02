clc
clear 
close all
warning off
n = 2; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 0.1;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
Sigma = fun_rho(rou,N,1,0.1);
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
for i =1:100
    Train = fun_TrainData(str_train,N,L,Sigma,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    [x0,tau] = fun_TrainData(str_train,N,1,Sigma,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_LogM = fun_RLogEMean(Train);
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [R_CC,alpha_cc(i)]=fun_CC(Train,R_NSCM,R_KA);
    [R_ML,alpha_ML(i)]=fun_MLalpha(Train,R_NSCM,R_KA,x0);
    [R_ECC,alpha_ecc(i)]=fun_ECC(Train,R_NSCM,R_KA);
    [R_LECC,alpha_lecc(i)]=fun_LECC(Train,R_LogM,R_KA);
    error_RCC(i) = norm(R_CC-Sigma,'fro');
    error_RECC(i) = norm(R_ECC-Sigma,'fro');
    error_RLECC(i) = norm(R_LECC-Sigma,'fro');
    error_RML(i) = norm(R_ML-Sigma,'fro');
    error_RSCM(i) = norm(R_SCM-Sigma,'fro');
    error_LogM(i) = norm(R_LogM-Sigma,'fro');
    error_RNSCM(i) = norm(R_NSCM-Sigma,'fro');
end
toc
m_errorRCC = mean(error_RCC)/norm(Sigma,'fro');
m_errorRECC = mean(error_RECC)/norm(Sigma,'fro');
m_errorRLECC = mean(error_RLECC)/norm(Sigma,'fro');
m_errorRML = mean(error_RML)/norm(Sigma,'fro');
m_errorRSCM = mean(error_RSCM)/norm(Sigma,'fro');
m_errorLogM = mean(error_LogM)/norm(Sigma,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(Sigma,'fro');
m_errorRKA = norm(R_KA-Sigma,'fro')/norm(Sigma,'fro');

mean_alpha_ecc = mean(alpha_ecc);
mean_alpha_lecc = mean(alpha_lecc);
mean_alpha_cc = mean(alpha_cc);
mean_alpha_ML = mean(alpha_ML);

% num = 1:1000;
% plot(num,alpha_cc,'b')
% hold on
% plot(num,alpha_ML,'g')
% plot(num,alpha_ecc,'r')
% plot(num,alpha_lecc,'k')
% legend('CC','ML','E','LE')




