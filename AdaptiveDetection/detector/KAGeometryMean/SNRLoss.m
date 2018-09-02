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
sigma_t =0.9;
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
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
rouR = fun_rho(rou,N,1);
rouR_half=rouR^0.5;
irouR=inv(rouR);
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
R_KA = zeros(size(rouR));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + rouR.*(t*t')/1000;
end


