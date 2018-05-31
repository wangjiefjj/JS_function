%% ͼ4.32���������ƫ���Ǻ�ƫ��������\alpha_1�����߹�ϵ(H=506km)
clc;clear;close all
%% ���ƽ̨��������1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [506,1000,2000,7000];%km
alpha1 = (0:90)/180*pi; %�������µ�����γ��
eta = 90/180*pi; % fai���״��������н�,90�ȡ������
%% ���㺽ƫ��, ��ƫ����
crabA1 = abs(fun_CrabAngle(alpha1,eta,H(1))/pi*180);
crabM1 = abs(fun_CrabMagnitude(alpha1,eta,H(1)));

crabA2 = abs(fun_CrabAngle(alpha1,eta,H(2))/pi*180);
crabM2 = abs(fun_CrabMagnitude(alpha1,eta,H(2)));

crabA3 = abs(fun_CrabAngle(alpha1,eta,H(3))/pi*180);
crabM3 = fun_CrabMagnitude(alpha1,eta,H(3));

crabA4 = abs(fun_CrabAngle(alpha1,eta,H(4))/pi*180);
crabM4 = fun_CrabMagnitude(alpha1,eta,H(4));
%% figure
figure(1)
hold on
plot(alpha1/pi*180,crabA1,'r','LineWidth',2);
plot(alpha1/pi*180,crabA2,'k','LineWidth',2);
plot(alpha1/pi*180,crabA3,'g','LineWidth',2);
plot(alpha1/pi*180,crabA4,'b','LineWidth',2);
h_leg = legend('H=506km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('γ�� \alpha_1/^o')
ylabel('��ƫ�� \phi_c/^o')
title('�������ƫ���Ǻ�ƫ��������\alpha_1�����߹�ϵ')
%% figure
figure(2)
hold on
plot(alpha1/pi*180,crabM1,'r','LineWidth',2);
plot(alpha1/pi*180,crabM2,'k','LineWidth',2);
plot(alpha1/pi*180,crabM3,'g','LineWidth',2);
plot(alpha1/pi*180,crabM4,'b','LineWidth',2);
h_leg = legend('H=506km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
grid on
box on
xlabel('γ�� \alpha_1/^o')
ylabel('��ƫ���� \rho_c/^o')
title('�������ƫ���Ǻ�ƫ��������\alpha_1�����߹�ϵ')