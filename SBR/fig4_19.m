%% 图4.19，主波束内距离模糊数和地距的关系
clc;clear;close all
%% 天基平台参数设置1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRF1 = 500; %Hz
Tr1 = 1/PRF1; %s
H = [500,1000,2000,7000]; %% 天基雷达距星下点的距离,km
EL_beam = 1/180*pi; % 俯仰向波束宽度。1deg->rad
AZ_beam = 1/180*pi; % 方位向波束宽度。1deg->rad
c = 3e5; %km 光速
%% 计算掠地角
[Graze1,R1,~] = fun_GrazeAngle(H(1));
[Graze2,R2,~] = fun_GrazeAngle(H(2));
[Graze3,R3,~] = fun_GrazeAngle(H(3));
[Graze4,R4,~] = fun_GrazeAngle(H(4));
Graze1=Graze1/180*pi;
Graze2=Graze2/180*pi;
Graze3=Graze3/180*pi;
Graze4=Graze4/180*pi;
%% 计算最大不模糊距离 式 4.35
dR11 =  c*Tr1/2 *sec(Graze1);
dR21 =  c*Tr1/2 *sec(Graze2);
dR31 =  c*Tr1/2 *sec(Graze3);
dR41 =  c*Tr1/2 *sec(Graze4);
%% 计算主足印长度
L1 = fun_FootLength(H(1),EL_beam);
L2 = fun_FootLength(H(2),EL_beam);
L3 = fun_FootLength(H(3),EL_beam);
L4 = fun_FootLength(H(4),EL_beam);
%% 计算主波束内的距离模糊数Na
Na11 = ceil(L1(500:end-500)./(dR11(500:end-500)));
Na21 = ceil(L2(500:end-600)./(dR21(500:end-600)));
Na31 = ceil(L3(500:end-750)./(dR31(500:end-750)));
Na41 = ceil(L4(500:end-1200)./(dR41(500:end-1200)));
%% figure
figure(1)
% suptitle('总标题')
set(gcf,'Position',[700 0 600 1000])
subplot(2,1,1)
hold on
plot(R1(500:end-500),Na11,'r','LineWidth',2);
plot(R2(500:end-600),Na21,'k','LineWidth',2);
plot(R3(500:end-750),Na31,'g','LineWidth',2);
plot(R4(500:end-1200),Na41,'b','LineWidth',2);
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
axis([500,3500,0,3])
grid on
box on
xlabel('地距R/km')
ylabel('Na')
title('图4.17（a）主波束呢距离模糊数Na与地距的关系, PRF=500Hz')
%% 天基平台参数设置2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRF2 = 2000; %Hz
Tr2 = 1/PRF2; %s
%% 计算最大不模糊距离 式 4.35
dR12 =  c*Tr2/2 *sec(Graze1);
dR22 =  c*Tr2/2 *sec(Graze2);
dR32 =  c*Tr2/2 *sec(Graze3);
dR42 =  c*Tr2/2 *sec(Graze4);
%% 计算主波束内的距离模糊数Na
Na12 = ceil(L1(500:end-500)./(dR12(500:end-500)));
Na22 = ceil(L2(500:end-600)./(dR22(500:end-600)));
Na32 = ceil(L3(500:end-750)./(dR32(500:end-750)));
Na42 = ceil(L4(500:end-1200)./(dR42(500:end-1200)));
%% figure
subplot(2,1,2)
hold on
plot(R1(500:end-500),Na12,'r','LineWidth',2);
plot(R2(500:end-600),Na22,'k','LineWidth',2);
plot(R3(500:end-750),Na32,'g','LineWidth',2);
plot(R4(500:end-1200),Na42,'b','LineWidth',2);
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
axis([500,3500,0,7])
grid on
box on
xlabel('地距R/km')
ylabel('Na')
title('图4.17（a）主波束呢距离模糊数Na与地距的关系, PRF=2000Hz')
