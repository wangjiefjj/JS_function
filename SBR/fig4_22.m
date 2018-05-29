%% 图4.22，总的模糊数和地距的关系
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRF = 0:2000; %Hz
Tr = 1./PRF; %s
H = [500,1000,2000,7000]; %% 天基雷达距星下点的距离,km
c = 3e5; %km 光速
%% 计算最大斜距,前向模糊到最大斜距处
Rsmax1 = fun_Rsmax(H(1));
Rsmax2 = fun_Rsmax(H(2));
Rsmax3 = fun_Rsmax(H(3));
Rsmax4 = fun_Rsmax(H(4));
%% 计算最大斜距对应的地距,前向模糊到最大斜距处
Rmax1 = fun_Rs2R(H(1),Rsmax1);
Rmax2 = fun_Rs2R(H(2),Rsmax2);
Rmax3 = fun_Rs2R(H(3),Rsmax3);
Rmax4 = fun_Rs2R(H(4),Rsmax4);
%% 假设目标斜距在Rsmax和H的中间 
Rst1 = (Rsmax1+H(1))/2;
Rst2 = (Rsmax2+H(2))/2;
Rst3 = (Rsmax3+H(3))/2;
Rst4 = (Rsmax4+H(4))/2;
%% 计算目标对应的地距 
Rt1 = fun_Rs2R(H(1),Rst1);
Rt2 = fun_Rs2R(H(2),Rst2);
Rt3 = fun_Rs2R(H(3),Rst3);
Rt4 = fun_Rs2R(H(4),Rst4);
%% 计算掠地角
Graze1 = fun_GrazeAngle(H(1),Rt1,Rst1)/180*pi;
Graze2 = fun_GrazeAngle(H(2),Rt2,Rst2)/180*pi;
Graze3 = fun_GrazeAngle(H(3),Rt3,Rst3)/180*pi;
Graze4 = fun_GrazeAngle(H(4),Rt4,Rst4)/180*pi;
%% 计算最大不模糊地距离 式 4.35
dR1 =  c*Tr/2 *sec(Graze1);
dR2 =  c*Tr/2 *sec(Graze2);
dR3 =  c*Tr/2 *sec(Graze3);
dR4 =  c*Tr/2 *sec(Graze4);
%% 计算总模糊数
Na1 = ceil(Rmax1./dR1);
Na2 = ceil(Rmax2./dR2);
Na3 = ceil(Rmax3./dR3);
Na4 = ceil(Rmax4./dR4);
%% figure
figure(1)
hold on
plot(PRF,Na1,'r','LineWidth',2);
plot(PRF,Na2,'k','LineWidth',2);
plot(PRF,Na3,'g','LineWidth',2);
plot(PRF,Na4,'b','LineWidth',2);
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
% axis([100,2000,0,50])
grid on
box on
xlabel('PRF/Hz')
ylabel('Na')
title('图4.22 距离模糊数Na与地距、天际雷达高度的关系')
