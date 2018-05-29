%% ͼ4.8����������ӡ���ȣ���������ӡ�����ؾ�Ĺ�ϵ
clc;clear;close all
%% ���ƽ̨����
H = [500,1000,2000,7000]; %% ����״�����µ�ľ���,km
EL_beam = 1/180*pi; % ����������ȡ�1deg->rad
AZ_beam = 1/180*pi; % ��λ������ȡ�1deg->rad
%%

[R1,Rs1] = fun_RsR(H(1));
[R2,Rs2] = fun_RsR(H(2));
[R3,Rs3] = fun_RsR(H(3));
[R4,Rs4] = fun_RsR(H(4));
%% ��������ӡ����
L1 = fun_FootLength(H(1),EL_beam);
L2 = fun_FootLength(H(2),EL_beam);
L3 = fun_FootLength(H(3),EL_beam);
L4 = fun_FootLength(H(4),EL_beam);
%% figure
figure(1)
hold on
plot(R1(1:end-500),L1(1:end-500),'r','LineWidth',2);
plot(R2(1:end-600),L2(1:end-600),'k','LineWidth',2);
plot(R3(1:end-750),L3(1:end-750),'g','LineWidth',2);
plot(R4(1:end-1200),L4(1:end-1200),'b','LineWidth',2);
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('�ؾ�R/km')
ylabel('��������ӡ����/km')
title('��������ӡ������ؾ�Ĺ�ϵ')
%% ������ӡ���
W1 = Rs1*AZ_beam;
W2 = Rs2*AZ_beam;
W3 = Rs3*AZ_beam;
W4 = Rs4*AZ_beam;
%% figure
figure(2)
hold on
plot(R1,W1,'r','LineWidth',2);
plot(R2,W2,'k','LineWidth',2);
plot(R3,W3,'g','LineWidth',2);
plot(R4,W4,'b','LineWidth',2);
h_leg = legend('H=500km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('�ؾ�R/km')
ylabel('��������ӡ���/km')
title('��������ӡ�����ؾ�Ĺ�ϵ')
