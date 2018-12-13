%%仿真数据估计逆gamma纹理参数
clc
clear 
close all
rho = 0.9;  %%协方差矩阵生成的迟滞因子
N=8;
L=1;
str_train = 'p';
opt_train = 1;     
%%杂波协方差
CNRout=30;
CNRnum=10.^(CNRout/10); 
fc=0.1;
Rc = fun_rho(rho,N,1,fc);
Rc = Rc+ 1/CNRnum * eye(N) ;%+ eye(N)
%%%
lambda=2;
mu=0.5;
for i = 1:1000%nsweep
    i
    r = fun_TrainData(str_train,1,L,1,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
    [lambda_t(i),mu_t(i)] = fun_IG_ML(r(1,:));
end
lambda_m=mean(lambda_t)
mu_m=mean(mu_t)