%% 图4.25，距离方位域等多普勒线
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 4.5e9; %载频Hz
c = 3e8;% m
lambda = c/fo;
PRF = 500;
Tr = 1/PRF;
H = 500;%km
angle_Az = (0:180)./180*pi; % 方位角
cosAz = cos(angle_Az);
R = 0:3000; %地距km
%% 计算多普勒
for i = 1:length(angle_Az)
    wd(i,:) = fun_Wd(H,R,angle_Az(i), Tr, lambda);
end
%% figure
figure(1)
contour(R,angle_Az/pi*180,wd,20)
xlabel('地距R/km')
ylabel('方位角/^{o}')
title('图4.25 距离方位域等多普勒线')
%% figure
figure(2)
contour(R,cosAz,wd,20)
xlabel('地距R/km')
ylabel('cos(\theta_{AZ})')
title('图4.25 距离方位域等多普勒线')
%% figure
figure(3)
[X,Y] = meshgrid(R,angle_Az/pi*180);
mesh(X,Y,wd)
xlabel('地距R/km')
ylabel('方位角/^{o}')
zlabel('多普勒/Hz')
title('图4.25 距离方位域等多普勒线')
