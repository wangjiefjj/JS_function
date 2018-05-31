%% 图4.5 最大地距，最大俯仰角与高度的关系
clc;clear;close all
%% 天基平台参数
H = 100:10000; %% 天基雷达距星下点的距离,km
%% 计算最大地距
Rmax = fun_Rmax(H);
%% figure
figure(1)
plot(H,Rmax ,'LineWidth',2)
grid on
box on
xlabel('高度H/km')
ylabel('最大地距Rmax/km')
title('图4.5 最大地距与高度的关系')
%% 计算最大俯仰角
ELmax = fun_ELmax(H);
%% figure
figure(2)
plot(H,ELmax ,'LineWidth',2)
grid on
box on
xlabel('高度H/km')
ylabel('最大俯仰角\theta_{EL,max}/km')
title('图4.5 最大俯仰角与高度的关系')
