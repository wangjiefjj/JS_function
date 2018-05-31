%% ͼ4.2��б����ؾ�Ĺ�ϵ
clc;clear;close all
%% ���ƽ̨����
H = [500,1000,2000,7000]; %% ����״�����µ�ľ���,km
%% ����ؾ�
[R1,Rs1] = fun_RsR(H(1));

[R2,Rs2] = fun_RsR(H(2));

[R3,Rs3] = fun_RsR(H(3));

[R4,Rs4] = fun_RsR(H(4));
%% Figure
hold on
plot(R1,Rs1,'r','LineWidth',2);
plot(R2,Rs2,'k','LineWidth',2);
plot(R3,Rs3,'g','LineWidth',2);
plot(R4,Rs4,'b','LineWidth',2);

plot(R1(end),Rs1(end),'ro')
plot(R1,Rs1(end)*ones(1,length(R1)),'k--','LineWidth',2)
plot(R1(end)*ones(1,length(R1)),linspace(0,Rs1(end),length(Rs1)),'k--','LineWidth',2)

plot(R2(end),Rs2(end),'ro')
plot(R2,Rs2(end)*ones(1,length(R2)),'k--','LineWidth',2)
plot(R2(end)*ones(1,length(R2)),linspace(0,Rs2(end),length(Rs2)),'k--','LineWidth',2)

plot(R3(end),Rs3(end),'ro')
plot(R3,Rs3(end)*ones(1,length(R3)),'k--','LineWidth',2)
plot(R3(end)*ones(1,length(R3)),linspace(0,Rs3(end),length(Rs3)),'k--','LineWidth',2)

plot(R4(end),Rs4(end),'ro')
plot(R4,Rs4(end)*ones(1,length(R4)),'k--','LineWidth',2)
plot(R4(end)*ones(1,length(R4)),linspace(0,Rs4(end),length(Rs4)),'k--','LineWidth',2)

grid on
box on
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
xlabel('�ؾ�R/km')
ylabel('б��R_s/km')
title('б����ؾ�Ĺ�ϵ')


