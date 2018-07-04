%% 图4.24，多普勒频率与距离和方位角的关系曲线
clc;clear;close all
%% 天基平台参数设置1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 4.5e9; %载频Hz
c = 3e8;% m
lambda = c/fo;
PRF = 500;
Tr = 1/PRF;
H = 500e3;%km
angle_Az = (0:180)./180*pi; % 方位角
cosAz = cos(angle_Az);
R = [250e3,1000e3,1750e3]; %地距km
%% 计算多普勒
wd1 = fun_Wd(H,R(1),angle_Az, Tr, lambda);
wd2 = fun_Wd(H,R(2),angle_Az, Tr, lambda);
wd3 = fun_Wd(H,R(3),angle_Az, Tr, lambda);
%% figure
figure(1)
set(gcf,'Position',[700 0 600 1000])
subplot(2,1,1)
hold on
plot(cosAz,wd1,'r','LineWidth',2);
plot(cosAz,wd2,'k','LineWidth',2);
plot(cosAz,wd3,'g','LineWidth',2);
h_leg = legend('R=250km','R=1000km','R=1750km');
set(h_leg,'Location','SouthEast')
% axis([100,2000,0,50])
grid on
box on
xlabel('cos(\theta_{AZ})')
ylabel('多普勒/Hz')
title('图4.24(a) 多普勒频率与距离和方位角的关系')
%% 天基平台参数设置2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
fo = 4.5e9; %载频Hz
c = 3e8;% m
lambda = c/fo;
PRF = 500;
Tr = 1/PRF;
H = 500;%km
angle_Az = [60,75,85,90]./180*pi; % 方位角
cosAz = cos(angle_Az);
R = 0:3000; %地距km
%% 计算多普勒
wd12 = fun_Wd(H,R,angle_Az(1), Tr, lambda);
wd22 = fun_Wd(H,R,angle_Az(2), Tr, lambda);
wd32 = fun_Wd(H,R,angle_Az(3), Tr, lambda);
wd42 = fun_Wd(H,R,angle_Az(4), Tr, lambda);
subplot(2,1,2)
hold on
plot(R,wd12,'r','LineWidth',2);
plot(R,wd22,'k','LineWidth',2);
plot(R,wd32,'g','LineWidth',2);
plot(R,wd42,'b','LineWidth',2);
h_leg = legend('\theta=60^{o}','\theta=75^{o}','\theta=85^{o}','\theta=90^{o}');
set(h_leg,'Location','SouthEast')
axis([0,2000,-5,120])
grid on
box on
xlabel('地距R/km')
ylabel('多普勒/Hz')
title('图4.24(b) 多普勒频率与距离和方位角的关系')
