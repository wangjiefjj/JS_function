%% ͼ4.3���ӵؽǣ���������ؾ�Ĺ�ϵ
clc;clear;close all
%% ���ƽ̨����
H = [500,1000,2000,7000]; %% ����״�����µ�ľ���,km
%% ���������
[GA1,R1,Rs1] = fun_GrazeAngle(H(1)); %%�ӵؽ�
[GA2,R2,Rs2] = fun_GrazeAngle(H(2)); %%�ӵؽ�
[GA3,R3,Rs3] = fun_GrazeAngle(H(3)); %%�ӵؽ�
[GA4,R4,Rs4] = fun_GrazeAngle(H(4)); %%�ӵؽ�

%% Figure
figure(1)
hold on
plot(R1,GA1,'r','LineWidth',2);
plot(R2,GA2,'k','LineWidth',2);
plot(R3,GA3,'g','LineWidth',2);
plot(R4,GA4,'b','LineWidth',2);
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('�ؾ�R/km')
ylabel('�ӵؽ�/deg')
title('�ӵؽ���ؾ�Ĺ�ϵ')
% %% ���㸩����
[EL1,R1,Rs1] = fun_ELAngle(H(1)); %%�ӵؽ�
[EL2,R2,Rs2] = fun_ELAngle(H(2)); %%�ӵؽ�
[EL3,R3,Rs3] = fun_ELAngle(H(3)); %%�ӵؽ�
[EL4,R4,Rs4] = fun_ELAngle(H(4)); %%�ӵؽ�
% %% Figure
figure(2)
hold on
plot(R1,EL1,'r','LineWidth',2);
plot(R2,EL2,'k','LineWidth',2);
plot(R3,EL3,'g','LineWidth',2);
plot(R4,EL4,'b','LineWidth',2);

plot(R1(end),EL1(end),'ro')
plot(R1,EL1(end)*ones(1,length(EL1)),'k--','LineWidth',2)
plot(R1(end)*ones(1,length(R1)),linspace(0,EL1(end),length(EL1)),'k--','LineWidth',2)

plot(R2(end),EL2(end),'ro')
plot(R2,EL2(end)*ones(1,length(EL2)),'k--','LineWidth',2)
plot(R2(end)*ones(1,length(R2)),linspace(0,EL2(end),length(EL2)),'k--','LineWidth',2)

plot(R3(end),EL3(end),'ro')
plot(R3,EL3(end)*ones(1,length(EL3)),'k--','LineWidth',2)
plot(R3(end)*ones(1,length(R3)),linspace(0,EL3(end),length(EL3)),'k--','LineWidth',2)

plot(R4(end),EL4(end),'ro')
plot(R4,EL4(end)*ones(1,length(EL4)),'k--','LineWidth',2)
plot(R4(end)*ones(1,length(R4)),linspace(0,EL4(end),length(EL4)),'k--','LineWidth',2)
grid on
box on
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
xlabel('�ؾ�R/km')
ylabel('������/deg')
title('��������ؾ�Ĺ�ϵ')
