clc
clear 
close all
% warning off
n = 0.5; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 2;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 0.1;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
CNRout=10; % 输出CNR
PFA=1e-3;% PFA=1e-4;
CNRnum=10.^(CNRout/10);
R = fun_rho(rou,N,1,0.1);
R = CNRnum*R + eye(N);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.5;
nn = 0:N-1;

iter = 10;
for i =1:1000
    warning off
    i;
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_SCM2 = R_SCM + eye(N);
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    R_MESS = fun_MESS(R_SCM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_RSCM(i) = norm(R_SCM-R,'fro');
    error_RSCM2(i) = norm(R_SCM2-R,'fro');
    error_RNSCM(i) = norm(R_NSCM-R,'fro');
    error_RMESS(i) = norm(R_MESS-R,'fro');   
end
m_errorRSCM = mean(error_RSCM)/norm(R,'fro');
m_errorRSCM2 = mean(error_RSCM2)/norm(R,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(R,'fro');
m_errorRMESS = mean(error_RMESS)/norm(R,'fro');
% rank(R_SCM)
% rank(R_MESS)