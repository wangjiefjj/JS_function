clc
clear 
close all
%%%%参数设置
n = 0.5; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
rou = 0.90;  %%协方差矩阵生成的迟滞因子
sigma_t =(0.9);
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 5;     % 脉冲数
N = Na*Np;
PFA=1e-3;% PFA=1e-4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
rouR = fun_rho(rou,N,1);
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
R_KA = zeros(size(rouR));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + rouR.*(t*t')/1000;
end
tic
R_SCM = zeros(N,N);
R_CC = zeros(N,N);
R_ML = zeros(N,N);
R_NSCM = zeros(N,N);
R_ECC = zeros(N,N);
R_LogM = zeros(N,N);
R_LogCC = zeros(N,N);
R_P = zeros(N,N);
R_PCC = zeros(N,N);
R_SFP = zeros(N,N);
MC = 100;
for i = 1:MC
    i
%     warning('off')
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
    %%协方差估计
    R_SCM = R_SCM + fun_SCMN(Train)/MC;
    R_CC = R_CC + fun_CC(Train,fun_NSCMN(Train),R_KA)/MC;
    R_ML = R_ML + fun_MLalpha(Train,fun_NSCMN(Train),R_KA,x0)/MC;
    R_NSCM = R_NSCM + fun_NSCMN(Train)/MC;
    R_ECC = R_ECC + fun_PowerCC(Train,R_KA,1,4)/MC;
    R_LogM = R_LogM + fun_RLogEMean(Train,4)/MC;
    R_LogCC = R_LogCC + fun_LogCC_new(Train,R_KA,4)/MC;
    R_P = R_P + fun_RPowerEMean(Train,-1,4)/MC;
    R_PCC = R_PCC + fun_PowerCC(Train,R_KA,-1,4)/MC;
    R_SFP = R_SFP + fun_SFP(Train,1)/MC;
end
[PSD_R,ft] = fun_PSD(rouR);
[PSD_SCM] = fun_PSD(R_SCM);
[PSD_CC] = fun_PSD(R_CC);
[PSD_ML] = fun_PSD(R_ML);
[PSD_NSCM] = fun_PSD(R_NSCM);
[PSD_ECC] = fun_PSD(R_ECC);
[PSD_LogM] = fun_PSD(R_LogM);
[PSD_LogCC] = fun_PSD(R_LogCC);
[PSD_P] = fun_PSD(R_P);
[PSD_PCC] = fun_PSD(R_PCC);
[PSD_SFP] = fun_PSD(R_SFP);
figure
plot(ft,circshift(10*log10(PSD_R/max(PSD_R)),10),'k');
hold on
plot(ft,circshift(10*log10(PSD_CC/max(PSD_R)),10),'c');
plot(ft,circshift(10*log10(PSD_ML/max(PSD_R)),10),'r');
plot(ft,circshift(10*log10(PSD_NSCM/max(PSD_R)),10),'g-');
plot(ft,circshift(10*log10(PSD_ECC/max(PSD_R)),10),'g-.');
plot(ft,circshift(10*log10(PSD_LogM/max(PSD_R)),10),'b-');
plot(ft,circshift(10*log10(PSD_LogCC/max(PSD_R)),10),'b-.');
plot(ft,circshift(10*log10(PSD_P/max(PSD_R)),10),'m-');
plot(ft,circshift(10*log10(PSD_PCC/max(PSD_R)),10),'m-.');
plot(ft,circshift(10*log10(PSD_SFP/max(PSD_R)),10),'y-.');
grid on
h_leg = legend('OPT','CC','ML','NSCM','ECC','LogM','LogCC','P','PCC','SFP');