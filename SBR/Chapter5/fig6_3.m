%% 图6.3，距离500km及距离混叠处的多普勒频率，有地球自转
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;                             %阵元个数
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 506;                            %SBR高度km
Vp = fun_Vp(H);                     %SBR速度km/s
d = 13.4;                           %归一化阵元间距
PRF = 2000;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
c = 3e5;                            %光速km/m
beta = 19.47;
beta0 = beta * d;
alpha1 = 30/180*pi;                  %卫星纬度
eta = 45/180*pi;                    %卫星倾角
%% 目标
R=500;                              %真实目标距离km
Rs = fun_R2Rs(H,R);
Graze = fun_GrazeAngle(H, R, Rs)/180*pi;
dR =  c*Tr/2 *sec(Graze);   
RR=R+(-2:4)*dR;                     %模糊距离 7个
%% 计算多普勒
c = -1:0.01:1;                       %入射锥角
Az = (0:180)/180*pi;
ELm = fun_ELAngle(H,RR)/180*pi;
% wd=0;
for i = 1:length(ELm)
    cosAz = c./sin(ELm(i));
    wd(i,:) = beta0*sin(ELm(i)).*cosAz;
end
CrabA = (fun_CrabAngle( alpha1,eta, H));
CrabM = fun_CrabMagnitude( alpha1,eta, H);
% wdr=0;
for i = 1:length(ELm)
    Az = (acos(c./sin(ELm(i))));
    wdr(i,:) = beta0*CrabM*sin(ELm(i)).*cos(Az+CrabA);
end
wdr = real(wdr)+imag(wdr);
%% figure
figure(1)
hold on
plot(c,wd','b')
plot(c,(wdr.'),'r')
