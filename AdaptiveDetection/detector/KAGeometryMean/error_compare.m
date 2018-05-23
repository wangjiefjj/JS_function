clc
clear 
close all
% warning off
n = 2; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = [0.01,0.1:0.1:10];
%  sigma_t = [0.01:0.01:1];
% sigma_t = [11:101];
L_s = length(sigma_t);
L_R = 100;
opt = 'k';
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
Sigma = fun_rho(rou,N,2);
L=round(n*N); 
m_errorRCC = zeros(1,L_s);
m_errorRCCIter = zeros(1,L_s);
m_errorRCCML= zeros(1,L_s);
m_errorR_2 = zeros(1,L_s);
m_errorRSCM = zeros(1,L_s);
m_errorRKA = zeros(1,L_s);
for i_s = 1:L_s
    i_s
    R_KA = zeros(size(Sigma));
    for i = 1:1000
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
        R_KA = R_KA + Sigma.*(t*t')/1000;
    end
    error_RCC = zeros(1,L_R);
    error_RCCIter = zeros(1,L_R);
    error_RML = zeros(1,L_R);
    error_R_2 = zeros(1,L_R);
    error_RSCM = zeros(1,L_R);
    parfor i =1:L_R
        Train = fun_TrainData(str_train,N,L,Sigma,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = fun_TrainData(str_train,N,1,Sigma,lambda,mu,opt_train); 
        R_SCM = (fun_SCMN(Train));
        R_NSCM = (fun_NSCMN(Train));
        R_LE = fun_RLogEMean(Train);
        R_x0 = (fun_SCMN(x0));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
        [R_CC,alpha_cc(i)]=fun_CC(Train,R_NSCM,R_KA);
        [R_ML,alpha_ML(i)]=fun_MLalpha(Train,R_NSCM,R_KA,x0);
        [R_ECC,alpha_ecc(i)]=fun_ECC(Train,R_NSCM,R_KA,2);
        [R_LECC,alpha_lecc(i)]=fun_test_LogEKA(Train,R_NSCM,R_KA);
% %         [R_LECC,alpha_lecc(i)]=fun_LECC(Train,R_NSCM,R_KA,2);
        error_RCC(i) = norm(R_CC-Sigma,'fro');
        error_RECC(i) = norm(R_ECC-Sigma,'fro');
        error_RLECC(i) = norm(R_LECC-Sigma,'fro');
        error_RML(i) = norm(R_ML-Sigma,'fro');
        error_RSCM(i) = norm(R_SCM-Sigma,'fro');
        error_RNSCM(i) = norm(R_NSCM-Sigma,'fro');
    end

    m_errorRCC(i_s) = mean(error_RCC)/norm(Sigma,'fro');
    m_errorRECC(i_s) = mean(error_RECC)/norm(Sigma,'fro');
    m_errorRLECC(i_s) = mean(error_RLECC)/norm(Sigma,'fro');
    m_errorRML(i_s) = mean(error_RML)/norm(Sigma,'fro');
    m_errorRSCM(i_s) = mean(error_RSCM)/norm(Sigma,'fro');
    m_errorRNSCM(i_s) = mean(error_RNSCM)/norm(Sigma,'fro');
    m_errorRKA(i_s) = norm(R_KA-Sigma,'fro')/norm(Sigma,'fro');

    mean_alpha_ecc(i_s) = mean(alpha_ecc);
    mean_alpha_lecc(i_s) = mean(alpha_lecc);
    mean_alpha_cc(i_s) = mean(alpha_cc);
    mean_alpha_ML(i_s) = mean(alpha_ML);
end

figure
hold on 
plot(sigma_t,m_errorRSCM,'r','LineWidth',2)%k--
plot(sigma_t,m_errorRCC,'b','LineWidth',2)%k-.
plot(sigma_t,m_errorRECC,'k','LineWidth',2)
plot(sigma_t,m_errorRLECC,'c','LineWidth',2)
plot(sigma_t,m_errorRML,'g','LineWidth',2)%k:
grid on
h_leg = legend('SCM','CC','KA-CE','KA-LCE','ML');
xlabel('\sigma^2','FontSize',10)
ylabel('Error','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')
% axis([0,10,0,100])
figure
hold on 
plot(sigma_t,mean_alpha_cc,'b','LineWidth',2)%k-.
plot(sigma_t,mean_alpha_ecc,'k','LineWidth',2)
plot(sigma_t,mean_alpha_lecc,'c','LineWidth',2)
plot(sigma_t,mean_alpha_ML,'g','LineWidth',2)%k:
grid on
h_leg = legend('CC','KA-CE','KA-LCE','ML');
xlabel('\sigma^2','FontSize',10)
ylabel('Error','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')