clc
clear 
close all
%%%%参数设置
n = 2; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
rou = 0.95;  %%协方差矩阵生成的迟滞因子
sigma_t =sqrt(0.9);
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=-5:1:25; % 输出SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rouR = zeros(N,N);  %%真实的杂波协方差
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
rouR = fun_rho(rou,N,1);
rouR_abs=abs(rouR);
rouR_half=rouR^0.5;
irouR=inv(rouR);
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
iter = 1000;
parfor i =1:iter
%     warning off
    i;
    Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); 
    %     %%先验协方差
    R_KA = zeros(N,N);
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA =  (rouR).*(t*t');
    %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCMN(Train));  
    R_CC = fun_CC(Train,R_SCM,R_KA);
    R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
    R_ECCT = fun_PowerCC(Train,R_KA,1,4);
    R_ECCS = fun_PowerCC(Train,R_KA,1,8);
    R_ECCP = fun_PowerCC(Train,R_KA,1,7);
    %%%%%%%%%%%%%%%%
    [SINR_SCM(:,i)]=fun_SINR(R_SCM,rouR);
    SINR_CC(:,i)=fun_SINR(R_CC,rouR);
    SINR_ML(:,i)=fun_SINR(R_ML,rouR);
    SINR_ECCT(:,i)=fun_SINR(R_ECCT,rouR);
    SINR_ECCS(:,i)=fun_SINR(R_ECCS,rouR);
    SINR_ECCP(:,i)=fun_SINR(R_ECCP,rouR);
end
% % NSCM_PSD=fun_PSD(mean_R_NSCM);
[R_PSD,ft]=fun_SINR(rouR,rouR);
mean_SINR_SCM=mean(SINR_SCM,2);
mean_SINR_CC=mean(SINR_CC,2);
mean_SINR_ML=mean(SINR_ML,2);
mean_SINR_ECCT=mean(SINR_ECCT,2);
mean_SINR_ECCS=mean(SINR_ECCS,2);
mean_SINR_ECCP=mean(SINR_ECCP,2);

figure()
hold on
plot(ft,R_PSD,'r');
plot(ft,mean_SINR_SCM,'g');
plot(ft,mean_SINR_CC,'b');
plot(ft,mean_SINR_ML,'c');
plot(ft,mean_SINR_ECCT,'k');
plot(ft,mean_SINR_ECCS,'k-*');
plot(ft,mean_SINR_ECCP,'k-o');
h_leg = legend('R','SCM','CC','ML','ECCT','ECCS','ECCP');
grid on
box on
str = [str_train,'_','SINR','_',num2str(n),'N','_s',num2str(sigma_t^2),'.mat'];
save (str); 