%% 图4.31，三类轨道下，偏航角，偏航幅度和距离的关系
clc;clear;close all
%% 天基平台参数设置1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = (0:4e4);%km
alpha11 = 0/180*pi; %卫星星下点所在纬度
eta1 = 90/180*pi; % fai：雷达轨道与赤道夹角
%% 计算航偏角, 航偏幅度
crabA1 = fun_CrabAngle(alpha11,eta1,H)/pi*180;
crabM1 = fun_CrabMagnitude(alpha11,eta1,H);
%% figure
figure(1)
suptitle('图4.31 三类轨道下，偏航角，偏航幅度和距离的关系')
set(gcf,'Position',[100 0 800 1000])
subplot(3,2,1)
plot(H,crabA1);
grid on
box on
xlabel('H/km')
ylabel('偏航角/^o')
legend('(a) \phi_c, 极轨道的SBR(\alpha_1=0^o, \eta_i=90^o)')
subplot(3,2,2)
plot(H,crabM1);
grid on
box on
xlabel('H/km')
ylabel('偏航幅度')
% axis([0,10000,1,1.8])
legend('(b) \rho_c, 赤道轨道的SBR(\alpha_1=0^o, \eta_i=90^o)')
%% 天基平台参数设置2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = (0:4e4);%km
alpha12 = 0/180*pi; %卫星星下点所在纬度
eta2 = 0/180*pi; % fai：雷达轨道与赤道夹角
%% 计算航偏角, 航偏幅度
crabA2 = fun_CrabAngle(alpha12,eta2,H)/pi*180;
crabM2 = fun_CrabMagnitude(alpha12,eta2,H);
%% figure
suptitle('图4.31 三类轨道下，偏航角，偏航幅度和距离的关系')
subplot(3,2,3)
plot(H,crabA2);
grid on
box on
xlabel('H/km')
ylabel('偏航角/^o')
legend('(c) \phi_c, 极轨道的SBR(\alpha_1=0^o, \eta_i=0^o)')
subplot(3,2,4)
plot(H,crabM2);
grid on
box on
xlabel('H/km')
ylabel('偏航幅度')
% axis([0,10000,1,1.8])
legend('(d) \rho_c, 倾斜轨道的SBR(\alpha_1=0^o, \eta_i=0^o)')
%% 天基平台参数设置3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = (0:4e4);%km
alpha13 = 0/180*pi; %卫星星下点所在纬度
eta3 = 45/180*pi; % fai：雷达轨道与赤道夹角
%% 计算航偏角, 航偏幅度
crabA3 = fun_CrabAngle(alpha13,eta3,H)/pi*180;
crabM3 = fun_CrabMagnitude(alpha13,eta3,H);
%% figure
suptitle('图4.31 三类轨道下，偏航角，偏航幅度和距离的关系')
subplot(3,2,5)
plot(H,crabA3);
grid on
box on
xlabel('H/km')
ylabel('偏航角/^o')
legend('(c) \phi_c, 极轨道的SBR(\alpha_1=0^o, \eta_i=45^o)')
subplot(3,2,6)
plot(H,crabM3);
grid on
box on
xlabel('H/km')
ylabel('偏航幅度')
% axis([0,10000,1,1.8])
legend('(d) \rho_c, 极轨道的SBR(\alpha_1=0^o, \eta_i=45^o)')