%% 最小可检测速度
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9600e6; %450e6   9600e6                 %载频 Hz
C = 299792458;                            %光速 m/s
lambda = C/fo;                       %波长 m
Nr = 32;                             %接收阵元个数
Nt = 1;                              %发射阵元个数
Np = 32;                             %相干脉冲数
D = Nt*Nr*Np;                          %自由度 
Re = 6373e3;                          %地球半径m
H = 700e3;  %700e3  %9e3                 %SBR高度 m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR速度m/s
end                
d = 1;%                             %归一化阵元间距
gamma = 16;                         %发射接收空间比率
PRF = 500;                         %脉冲重复频率Hz，500~2000
Tr =1/PRF;                         %脉冲重复间隔s
% beta = 19.47;
beta = Vp*2/PRF/(d*lambda/2);      %%混叠系数 
alpha1 = 20/180*pi;                %卫星纬度
eta = 90/180*pi;                   %卫星倾角
%% 距离
Rs = 130d3;
if H>=500e3
    R = 1301e3;                               %地距m
    Rs = fun_R2Rs(H, R);
end

%% 距离模糊
% deltaR = fun_deltaR(H, R, Rs, Tr);
% R = [R];  %,R+2*dR,R+3*dR,R+4*dR
%%
EL = asin(H/Rs);         
if H>=500e3
    EL_du = fun_ELAngle(H,R);             %波束俯仰角
    EL = EL_du/180*pi;
end

%% 偏航角幅度
CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度
%% 杂波协方差
Num = 360;
CNR = 40;
Rk = zeros(D,D);
Rk_r = zeros(D,D);
for i = 1:length(EL)
    Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,lambda,gamma);%,CrabA,CrabM
    Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma,lambda,CrabA,CrabM);    
end
figure()
mesh(abs(Rk))
iRk = inv(Rk);  
iRk_r = inv(Rk_r); 

Train = fun_TrainData('p',D,2*D,Rk,2,0.5,1);
RSCM = fun_NSCMN(Train);
iRSCM = inv(RSCM);

%% 目标
wd = linspace(-0.5,0.5,200);               %虚拟规1化多普勒
wd_real = wd*PRF/2*lambda/2;                %虚拟实际速度m/s
cmj = wd./beta/(d*lambda^2);
fsp = d*cmj;

for i_wd = 1:length(wd)
    i_wd
    a = exp(1j*(0:Nr-1)*2*pi*fsp(i_wd)).';               %接收空间导向矢量
    b = exp(1j*(0:Nt-1)*2*pi*gamma*fsp(i_wd)).';         %发射空间导向矢量
    c = exp(1j*(0:Np-1)*2*pi*wd(i_wd)).';            %时间导向矢量
    s = kron(c,kron(b,a));                                  %遍历的导向矢量
    SINR_ideal(i_wd) = (abs(s'*iRk*s))^2;          % SINR_ideal 式6.37   
    SINR_ideal_r(i_wd) = (abs(s'*iRk_r*s))^2;      % SINR_ideal 式6.37   
    SINR(i_wd) = (abs(s'*iRSCM*s))^2/abs(s'*iRSCM*Rk*iRSCM*s);
end
% figure()
% plot(wd,10*log10(SINR_ideal/max(SINR_ideal)))%/max(SINR_ideal)
% hold on
% plot(wd,10*log10(SINR_ideal_r/max(SINR_ideal_r)))%/max(SINR_ideal)
% plot(wd,10*log10(SINR/max(SINR)),'r')%/max(SINR)

