clc
clear 
close all
% warning off
n = 1; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = sqrt([0.01,0.1,0.5,0.9]);
L_s = length(sigma_t);
L_R = 1000;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
R = fun_rho(rou,N,2);
L=round(n*N); 
m_errorRCC = zeros(1,L_s);
m_errorRECCT = zeros(1,L_s);
m_errorRECCS = zeros(1,L_s);
m_errorRECCP = zeros(1,L_s);
m_errorRML= zeros(1,L_s);
m_errorRSCM = zeros(1,L_s);
for i_s = 1:L_s
    i_s
    error_RCC = zeros(1,L_R);
    error_RECCT = zeros(1,L_R);
    error_RML = zeros(1,L_R);
    error_R_2 = zeros(1,L_R);
    error_RSCM = zeros(1,L_R);
    parfor i =1:L_R
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
        R_real = (tau0^2)*R;
        %%先验协方差
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
        R_KA =  (R).*(t*t');
        %%%
%         Train = Train/(diag(sqrt(diag(Train'*Train))));%%%归一化 
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCMN(Train));
        R_CC = fun_CC(Train,R_SCM,R_KA);
        R_ECCT = fun_PowerCC(Train,R_KA,1,4);
        R_ECCS = fun_PowerCC(Train,R_KA,1,8);
        R_ECCP = fun_PowerCC(Train,R_KA,1,7);
        R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
        %%NFN
        error_RCC(i) = norm(R_CC-R,'fro')/norm(R,'fro');
        error_RECCT(i) = norm(R_ECCT-R,'fro')/norm(R,'fro');
        error_RECCS(i) = norm(R_ECCS-R,'fro')/norm(R,'fro');
        error_RECCP(i) = norm(R_ECCP-R,'fro')/norm(R,'fro');
        error_RML(i) = norm(R_ML-R,'fro')/norm(R,'fro');    
        error_RSCM(i) = norm(R_SCM-R,'fro')/norm(R,'fro');
    end

    m_errorRCC(i_s) = mean(error_RCC);
    m_errorRECCT(i_s) = mean(error_RECCT);
    m_errorRECCS(i_s) = mean(error_RECCS);
    m_errorRECCP(i_s) = mean(error_RECCP);
    m_errorRML(i_s) = mean(error_RML);
    m_errorRSCM(i_s) = mean(error_RSCM);
end

figure
hold on 
plot(sigma_t,m_errorRSCM,'g','linewidth',2)
plot(sigma_t,m_errorRCC,'b','linewidth',2)
plot(sigma_t,m_errorRML,'c','linewidth',2)
plot(sigma_t,m_errorRECCT,'k','linewidth',2)
plot(sigma_t,m_errorRECCS,'K-*','linewidth',2)
plot(sigma_t,m_errorRECCP,'k-o','linewidth',2)
grid on
h_leg = legend('ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with T-KA','ANMF with S-KA','ANMF with P-KA');
xlabel('\sigma^2','FontSize',10)
ylabel('Error','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')
str = [str_train,'_Rerror_',num2str(n),'N','.mat'];
save (str); 