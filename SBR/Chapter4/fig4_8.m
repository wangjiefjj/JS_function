%% ͼ4.8����������ӡ���ȣ���������ӡ�����ؾ�Ĺ�ϵ
clc;clear;close all
%% ���ƽ̨����
H = [700,1000,2000,7000]; %% ����״�����µ�ľ���,km
EL_beam = 1.15/180*pi; % ����������ȡ�1deg->rad
AZ_beam = 0.57/180*pi; % ��λ������ȡ�1deg->rad
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
plot(R1(1:end-600),L1(1:end-600),'r','LineWidth',2);
plot(R2(1:end-700),L2(1:end-700),'k','LineWidth',2);
plot(R3(1:end-800),L3(1:end-800),'g','LineWidth',2);
plot(R4(1:end-1200),L4(1:end-1200),'b','LineWidth',2);
h_leg = legend('H=700km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('�ؾ�R/km')
ylabel('��������ӡ����/km')
title(['��������ӡ������ؾ�Ĺ�ϵ, \phi_{EL}=',num2str(EL_beam/pi*180),...
    '^{o}, \phi_{AZ}=',num2str(AZ_beam/pi*180),'^{o}'])
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
h_leg = legend('H=700km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('�ؾ�R/km')
ylabel('��������ӡ���/km')
title(['��������ӡ�����ؾ�Ĺ�ϵ, \phi_{EL}=',num2str(EL_beam/pi*180),...
    '^{o}, \phi_{AZ}=',num2str(AZ_beam/pi*180),'^{o}'])
