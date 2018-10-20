clc
clear
close all
clc
clear 
close all
% warning off
n = 1; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 2;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = 0.1;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
CNRout=30; % 输出CNR
PFA=1e-3;% PFA=1e-4;
CNRnum=10.^(CNRout/10);
R = fun_rho(rou,N,1,0);
R = CNRnum*R+eye(N);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.5;
nn = 0:N-1;
iter = 100;
mean_R_SCM=0;
mean_R_SCM2=0;
mean_R_MESS=0;
mean_R_NSCM=0;
for i =1:iter
%     warning off
    i;
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_SCM2 = R_SCM + eye(N);
    R_NSCM = (fun_NSCMN(Train));
    R_SFP = fun_SFP(Train,1);
    R_MESS = fun_MESS(R_SCM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_RSCM(i) = norm(R_SCM-R,'fro');
    error_RSCM2(i) = norm(R_SCM2-R,'fro');
    error_RNSCM(i) = norm(R_NSCM-R,'fro');
    error_RMESS(i) = norm(R_MESS-R,'fro');  
    %%%%%%%%%%%%%%%%
    [SCM_PSD(:,i)]=fun_SINR(R_SCM,R);
    SCM2_PSD(:,i)=fun_SINR(R_SCM2,R);
    MESS_PSD(:,i)=fun_SINR(R_MESS,R);
    NSCM_PSD(:,i)=fun_SINR(R_NSCM,R);
    SFP_PSD(:,i)=fun_SINR(R_SFP,R);
end
m_errorRSCM = mean(error_RSCM)/norm(R,'fro');
m_errorRSCM2 = mean(error_RSCM2)/norm(R,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(R,'fro');
m_errorRMESS = mean(error_RMESS)/norm(R,'fro');
% [R_PSD,ft]=fun_PSD(R);
% % SCM_PSD=fun_PSD(mean_R_SCM);
% % SCM2_PSD=fun_PSD(mean_R_SCM2);
% % MESS_PSD=fun_PSD(mean_R_MESS);
% % NSCM_PSD=fun_PSD(mean_R_NSCM);
[R_PSD,ft]=fun_SINR(R,R);
mean_SCM_PSD=mean(SCM_PSD,2);
mean_SCM2_PSD=mean(SCM2_PSD,2);
mean_NSCM_PSD=mean(NSCM_PSD,2);
mean_MESS_PSD=mean(MESS_PSD,2);
mean_SFP_PSD=mean(SFP_PSD,2);

figure()
hold on
plot(ft,R_PSD,'r');
plot(ft,mean_SCM_PSD,'g');
plot(ft,mean_SCM2_PSD,'b');
plot(ft,mean_NSCM_PSD,'k');
plot(ft,mean_MESS_PSD,'y');
plot(ft,mean_SFP_PSD,'c');
grid on
box on