clc
clear 
close all
% warning off
n = 1; %几倍的样本
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
R = fun_rho(rou,N,2);
L=round(n*N); 
m_errorRCC = zeros(1,L_s);
m_errorRCCIter = zeros(1,L_s);
m_errorRCCML= zeros(1,L_s);
m_errorR_2 = zeros(1,L_s);
m_errorRSCM = zeros(1,L_s);
m_errorRKA = zeros(1,L_s);
for i_s = 1:L_s
    i_s
%     for i = 1:1000
%         t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
%         R_KA = R_KA + R.*(t*t')/1000;
%     end
    error_RCC = zeros(1,L_R);
    error_RCCIter = zeros(1,L_R);
    error_RCCML = zeros(1,L_R);
    error_R_2 = zeros(1,L_R);
    error_RSCM = zeros(1,L_R);
    
    parfor i =1:L_R
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
        R_real = (tau0^2)*R;
        %%先验协方差
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
        R_KA =  (tau0^2*R).*(t*t');
        %%%
        R_SCM = (fun_SCMN(Train));
        R_LogM = fun_RLogEMean(Train);
        R_NSCM = (fun_NSCMN(Train));
        R_x0 = (fun_SCMN(x0));
        R_AML = fun_AML(Train);
        [R_CC,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
        [R_CCML,alpha_ML(i)]=fun_MLalpha(Train,R_SCM,R_KA,x0);
        [R_CCIter,alpha_iter(i)]=fun_CCIter2(Train,R_SCM,R_KA);
        [R_LECC, alpha_lecc(i)] = fun_LECC(Train,R_LogM,R_KA,2);
        
        error_RCC(i) = norm(R_CC-R_real,'fro')/norm(R_real,'fro');
        error_RCCIter(i) = norm(R_CCIter-R_real,'fro')/norm(R_real,'fro');
        error_RLECC(i) = norm(R_LECC-R_real,'fro')/norm(R_real,'fro');
        error_RCCML(i) = norm(R_CCML-R_real,'fro')/norm(R_real,'fro');    
        error_RSCM(i) = norm(R_SCM-R_real,'fro')/norm(R_real,'fro');
    end

    m_errorRCC(i_s) = mean(error_RCC);
    m_errorRCCIter(i_s) = mean(error_RCCIter);
    m_errorRLECC(i_s) = mean(error_RLECC);
    m_errorRCCML(i_s) = mean(error_RCCML);
    m_errorRSCM(i_s) = mean(error_RSCM);
    

    m_alpha(i_s) = mean(alpha);
    m_alpha_iter(i_s) = mean(alpha_iter);
    m_alpha_lecc(i_s) = mean(alpha_lecc);
    m_alpha_ML(i_s) = mean(alpha_ML);
end

figure
hold on 
plot(sigma_t,m_errorRSCM,'r','LineWidth',2)%k--
plot(sigma_t,m_errorRCC,'b','LineWidth',2)%k-.
plot(sigma_t,m_errorRCCIter, 'k','LineWidth',2)
plot(sigma_t,m_errorRCCML,'g','LineWidth',2)%k:
plot(sigma_t,m_errorRLECC,'y','LineWidth',2)%k:
grid on
% h_leg = legend('SCM','CC','KA-CE','ML');
xlabel('\sigma^2','FontSize',10)
ylabel('Error','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')
% axis([0,10,0,100])
