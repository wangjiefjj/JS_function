clc
clear 
close all
% warning off
n = 1; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 0.99;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
R = fun_rho(rou,N,1,0.00);
SNRout=-5:1:20; % 输出SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
R_KA = zeros(size(R));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + R.*(t*t')/1000;
end

iter = 10;
for i =1:1000
    warning off
    i
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
%     R_KA = trace(R_SCM)/N*eye(N);
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    R_LogMean = fun_RLogEMean(Train,2);
    R_PowerMean = fun_RPowerEMean(Train,2);%-1,0.5,2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [R_LCC,alcc(i)] = fun_LogCC_new(Train,R_KA);
    [R_PCC,apcc(i)] = fun_PowerCC(Train,R_KA,2);%-1,0.5,2
    [R_CC,alpha(i)]=fun_CC(Train,R_NSCM,R_KA);
    error_RKA(i) = norm(R_KA-R,'fro');
    error_RLCC(i) = norm(R_LCC-R,'fro');
    error_RPCC(i) = norm(R_PCC-R,'fro');
    error_RP(i) = norm(R_PowerMean-R,'fro');
    error_RL(i) = norm(R_LogMean-R,'fro');
    error_RCC(i) = norm(R_CC-R,'fro');
    error_RSCM(i) = norm(R_SCM-R,'fro');
    error_RNSCM(i) = norm(R_NSCM-R,'fro');
    
end
m_errorRLCC = mean(error_RLCC)/norm(R,'fro');
m_errorRPCC = mean(error_RPCC)/norm(R,'fro');
m_errorRP = mean(error_RP)/norm(R,'fro');
m_errorRL = mean(error_RL)/norm(R,'fro');
m_errorRCC = mean(error_RCC)/norm(R,'fro');
m_errorRSCM = mean(error_RSCM)/norm(R,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(R,'fro');
m_errorRKA = mean(error_RKA)/norm(R,'fro');

mean_alcc = mean(alcc);
mean_apcc = mean(apcc);
mean_alpha = mean(alpha);

