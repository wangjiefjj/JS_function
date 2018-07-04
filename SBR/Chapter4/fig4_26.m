%% 图4.25，距离及杂波斜率(\beta_0=6)引起的多普勒扩展
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1.25e9; %载频Hz
c = 3e8;% m
lambda = c/fo;
H = 500e3;%m
Vp = fun_Vp(H);
angle_Az = (0:180)./180*pi; % 方位角
R = [1000e3,500e3]; %地距km
beta0 = 6;
Tr =  lambda/2*beta0/2/Vp;
wd1 = mod(fun_Wd_beta0(H,R(1),angle_Az,beta0),2)-1;
wd2 = mod(fun_Wd_beta0(H,R(2),angle_Az,beta0),2)-1;
%% figure
figure(1)
hold on
plot(cos(angle_Az),wd1,'r','LineWidth',2);
plot(cos(angle_Az),wd2,'k','LineWidth',2);
h_leg = legend('R=1000km','R=500km');
set(h_leg,'Location','SouthEast')
% axis([100,2000,0,50])
grid on
box on
xlabel('Azimuth')
ylabel('多普勒/Hz')
title('图4.26 距离及杂波斜率(\beta_0=6)引起的多普勒扩展')