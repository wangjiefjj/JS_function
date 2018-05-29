%% 图4.33，倾斜轨道，偏航角和偏航幅度与纬度的关系(H=506km)
clc;clear;close all
%% 天基平台参数设置1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [506,1000,7000];
alpha1 = (0:80)/180*pi; %卫星星下点所在纬度
eta = (0:80)/180*pi; % fai：雷达轨道与赤道夹角
[X,Y] = meshgrid(eta,alpha1);
%% 计算航偏角, 航偏幅度
crabA1 = real(fun_CrabAngle(X,Y,H(1))/pi*180);
crabM1 = (fun_CrabMagnitude(X,Y,H(1)));
%% figure
figure(1)
mesh(X/pi*180,Y/pi*180,crabA1);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\phi_c')
title('倾斜轨道，偏航角与纬度的关系(H=506km)')
%% figure
figure(2)
mesh(X/pi*180,Y/pi*180,crabM1);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\rho_c')
title('倾斜轨道，偏航幅度与纬度的关系(H=506km)')

%% 计算航偏角, 航偏幅度
crabA2 = real(fun_CrabAngle(X,Y,H(2))/pi*180);
crabM2 = (fun_CrabMagnitude(X,Y,H(2)));
%% figure
figure(3)
mesh(X/pi*180,Y/pi*180,crabA2);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\phi_c')
title('倾斜轨道，偏航角与纬度的关系(H=1000km)')
%% figure
figure(4)
mesh(X/pi*180,Y/pi*180,crabM2);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\rho_c')
title('倾斜轨道，偏航幅度与纬度的关系(H=1000km)')


%% 计算航偏角, 航偏幅度
crabA3 = real(fun_CrabAngle(X,Y,H(3))/pi*180);
crabM3 = (fun_CrabMagnitude(X,Y,H(3)));
%% figure
figure(5)
mesh(X/pi*180,Y/pi*180,crabA3);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\phi_c')
title('倾斜轨道，偏航角与纬度的关系(H=7000km)')
%% figure
figure(6)
mesh(X/pi*180,Y/pi*180,crabM3);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\rho_c')
title('倾斜轨道，偏航幅度与纬度的关系(H=7000km)')