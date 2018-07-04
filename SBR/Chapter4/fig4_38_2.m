%% 图4.38，偏航角对多普勒频率与距离的影响,mesh 图，方位角-30:30
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 12.5e9;                                  %载频 Hz
c = 3e8;                                    %光速 m/s
lambda = c/fo;                              %波长 m
PRF = 500;                                  %脉冲重复频率 Hz
Tr = 1/PRF;                                 %脉冲重复间隔 s
H = 500;                                    %卫星高度 km
angle_Az = [-90:1:90]./180*pi;              %方位角 rad
alpha1 = 30/180*pi;                         %纬度 rad  
eta = 45/180*pi;                            %轨道倾角 rad
R = 0:3000;                                 %地距 km
for i = 1:length(angle_Az)
%% 计算无自转多普勒
    wd(i,:) = fun_Wd(H,R,angle_Az(i), Tr, lambda);
    %% 计算有自转多普勒
    wdr(i,:) = fun_Wd_Rotation(H,R,angle_Az(i), Tr, lambda, alpha1, eta);
end
%% figure
figure(1)
[X,Y] = meshgrid(R,angle_Az./pi.*180);
mesh(X,Y, abs(wdr))
hold on
mesh(X,Y, abs(wd))
% axis([0,3000,-15,60])
grid on
box on
xlabel('地距/km')
ylabel('方位角/deg')
zlabel('多普勒/Hz')
title(['偏航角对多普勒频率与距离的影响(H=',num2str(H),'km)'])
