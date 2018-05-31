%% 图4.38，偏航角对多普勒频率与距离的影响
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 40e9;                                  %载频 Hz
c = 3e8;                                    %光速 m/s
lambda = c/fo;                              %波长 m
PRF = 500;                                  %脉冲重复频率 Hz
Tr = 1/PRF;                                 %脉冲重复间隔 s
H = 500;                                    %卫星高度 km
angle_Az = [75,80,90]./180*pi;              %方位角 rad
alpha1 = 30/180*pi;                         %纬度 rad  
eta = 45/180*pi;                            %轨道倾角 rad
R = 0:3000;                                 %地距 km
%% 计算无自转多普勒
wd1 = fun_Wd(H,R,angle_Az(1), Tr, lambda);
wd2 = fun_Wd(H,R,angle_Az(2), Tr, lambda);
wd3 = fun_Wd(H,R,angle_Az(3), Tr, lambda);
%% 计算有自转多普勒
wdr1 = fun_Wd_Rotation(H,R,angle_Az(1), Tr, lambda, alpha1, eta);
wdr3 = fun_Wd_Rotation(H,R,angle_Az(2), Tr, lambda, alpha1, eta);
wdr4 = fun_Wd_Rotation(H,R,angle_Az(3), Tr, lambda, alpha1, eta);
%% figure
figure(1)
hold on

plot(R,wd1,'r--','LineWidth',2);
plot(R,wdr1,'r','LineWidth',2);

plot(R,wd2,'g--','LineWidth',2);
plot(R,wdr3,'g','LineWidth',2);

plot(R,wd3,'b--','LineWidth',2);
plot(R,wdr4,'b','LineWidth',2);

h_leg = legend('无自转\theta_{AZ}=75^{o}','考虑自转\theta_{AZ}=75^{o}',...
    '无自转\theta_{AZ}=80^{o}','考虑自转\theta_{AZ}=80^{o}',...
    '无自转\theta_{AZ}=90^{o}','考虑自转\theta_{AZ}=90^{o}');
set(h_leg,'Location','SouthEast')
axis([0,3000,-15,60])
grid on
box on
xlabel('地距R/km')
ylabel('多普勒/Hz')
title('图4.38 偏航角对多普勒频率与距离的影响(H=500km)')
