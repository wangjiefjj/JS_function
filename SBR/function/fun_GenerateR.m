function [ Rk ] = fun_GenerateR( Nt, Nr, Np, Num, CNR, EL, beta, d,gamma, CrabA, CrabM)
%FUN_GENERATER 此处显示有关此函数的摘要
%   此处显示详细说明
%% 产生协方差矩阵, 假设无距离模糊, 可有地球自转
% Nt: 发射天线数
% Nr： 接收天线数
%Np： 脉冲数
% Num： 杂波块分数
% CNR: 杂噪比 
% EL: 俯仰角，rad
% CrabA: 航偏角
% CrabM：航片幅度
if nargin<10    %默认地球无自转
   CrabA = 0;
   CrabM = 1;
end
D = Nt*Nr*Np;
Rk = zeros(D,D);
AZ_du = linspace(0,179,Num);               %波束杂波块
% AZ_du = 90;
AZ = AZ_du/180*pi;
%% 阵列因子
steering_angle = 0; %% 波束指向
win = chebwin(Num);
for k=1:length(AZ)  
    AF(k) = sum(exp(-1i*pi*d*(0:Nr*Nt-1)*(sin(AZ(k)) ...
                  - sin(steering_angle*pi/180))));
end
AF(1) = AF(2);

AF = diag(abs(AF));

cmj = sin(EL)*cos(AZ);       %入射锥角
fspc = 0.5*d*cmj;             %空间频率
wdc = 0.5*beta*d * CrabM * sin(EL)*cos(AZ+CrabA);   %杂波归1化多普勒
sc = zeros(D,length(cmj));
for i_Az = 1:length(cmj)                   %俯仰角          
    a = exp(1j*(0:Nr-1)*2*pi*fspc(i_Az)).';         %接收空间导向矢量 
    b = exp(1j*(0:Nt-1)*2*pi*gamma*fspc(i_Az)).';    %发射空间导向矢量  
    c = exp(1j*(0:Np-1)*2*pi*wdc(i_Az)).';           %时间导向矢量 
    sc(:,i_Az) =  kron(c,kron(b,a));             %杂波的导向矢量%AF(i_Az) *
end 
Rk = sc*AF*sc';%*AF
Ac=(10^(CNR/10))^0.5; % 设噪声功率为1
Rk=D*Rk/sum(eig(Rk))*Ac^2;
Rk=Rk+eye(D); 
end

