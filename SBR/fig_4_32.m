%% 图4.32，极轨道，偏航角和偏航幅度与\alpha_1的曲线关系(H=506km)
clc;clear;close all
%% 天基平台参数设置1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [506,1000,2000,7000];%km
alpha1 = (0:90)/180*pi; %卫星星下点所在纬度
eta = 90/180*pi; % fai：雷达轨道与赤道夹角,90度·极轨道
%% 计算航偏角, 航偏幅度
crabA1 = abs(fun_CrabAngle(alpha1,eta,H(1))/pi*180);
crabM1 = abs(fun_CrabMagnitude(alpha1,eta,H(1)));

crabA2 = abs(fun_CrabAngle(alpha1,eta,H(2))/pi*180);
crabM2 = abs(fun_CrabMagnitude(alpha1,eta,H(2)));

crabA3 = abs(fun_CrabAngle(alpha1,eta,H(3))/pi*180);
crabM3 = fun_CrabMagnitude(alpha1,eta,H(3));

crabA4 = abs(fun_CrabAngle(alpha1,eta,H(4))/pi*180);
crabM4 = fun_CrabMagnitude(alpha1,eta,H(4));
%% figure
figure(1)
hold on
plot(alpha1/pi*180,crabA1,'r','LineWidth',2);
plot(alpha1/pi*180,crabA2,'k','LineWidth',2);
plot(alpha1/pi*180,crabA3,'g','LineWidth',2);
plot(alpha1/pi*180,crabA4,'b','LineWidth',2);
h_leg = legend('H=506km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('纬度 \alpha_1/^o')
ylabel('航偏角 \phi_c/^o')
title('极轨道，偏航角和偏航幅度与\alpha_1的曲线关系')
%% figure
figure(2)
hold on
plot(alpha1/pi*180,crabM1,'r','LineWidth',2);
plot(alpha1/pi*180,crabM2,'k','LineWidth',2);
plot(alpha1/pi*180,crabM3,'g','LineWidth',2);
plot(alpha1/pi*180,crabM4,'b','LineWidth',2);
h_leg = legend('H=506km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('纬度 \alpha_1/^o')
ylabel('航偏幅度 \rho_c/^o')
title('极轨道，偏航角和偏航幅度与\alpha_1的曲线关系')