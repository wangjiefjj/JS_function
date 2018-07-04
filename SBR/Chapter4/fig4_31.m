%% ͼ4.31���������£�ƫ���ǣ�ƫ�����Ⱥ;���Ĺ�ϵ
clc;clear;close all
%% ���ƽ̨��������1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = (0:10e3:4000e3);%m
alpha11 = 0/180*pi; %�������µ�����γ��
eta1 = 90/180*pi; % fai���״��������н�
%% ���㺽ƫ��, ��ƫ����
crabA1 = fun_CrabAngle(alpha11,eta1,H)/pi*180;
crabM1 = fun_CrabMagnitude(alpha11,eta1,H);
%% figure
figure(1)
suptitle('ͼ4.31 �������£�ƫ���ǣ�ƫ�����Ⱥ;���Ĺ�ϵ')
set(gcf,'Position',[100 0 800 1000])
subplot(3,2,1)
plot(H/1e3,crabA1);
grid on
box on
xlabel('H/km')
ylabel('ƫ����/^o')
legend('(a) \phi_c, �������SBR(\alpha_1=0^o, \eta_i=90^o)')
subplot(3,2,2)
plot(H/1e3,crabM1);
grid on
box on
xlabel('H/km')
ylabel('ƫ������')
% axis([0,10000,1,1.8])
legend('(b) \rho_c, ��������SBR(\alpha_1=0^o, \eta_i=90^o)')
%% ���ƽ̨��������2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = (0:10e3:4000e3);%m
alpha12 = 0/180*pi; %�������µ�����γ��
eta2 = 0/180*pi; % fai���״��������н�
%% ���㺽ƫ��, ��ƫ����
crabA2 = fun_CrabAngle(alpha12,eta2,H)/pi*180;
crabM2 = fun_CrabMagnitude(alpha12,eta2,H);
%% figure
suptitle('ͼ4.31 �������£�ƫ���ǣ�ƫ�����Ⱥ;���Ĺ�ϵ')
subplot(3,2,3)
plot(H/1e3,crabA2);
grid on
box on
xlabel('H/km')
ylabel('ƫ����/^o')
legend('(c) \phi_c, �������SBR(\alpha_1=0^o, \eta_i=0^o)')
subplot(3,2,4)
plot(H/1e3,crabM2);
grid on
box on
xlabel('H/km')
ylabel('ƫ������')
% axis([0,10000,1,1.8])
legend('(d) \rho_c, ��б�����SBR(\alpha_1=0^o, \eta_i=0^o)')
%% ���ƽ̨��������3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = (0:10e3:4000e3);%m
alpha13 = 0/180*pi; %�������µ�����γ��
eta3 = 45/180*pi; % fai���״��������н�
%% ���㺽ƫ��, ��ƫ����
crabA3 = fun_CrabAngle(alpha13,eta3,H)/pi*180;
crabM3 = fun_CrabMagnitude(alpha13,eta3,H);
%% figure
suptitle('ͼ4.31 �������£�ƫ���ǣ�ƫ�����Ⱥ;���Ĺ�ϵ')
subplot(3,2,5)
plot(H/1e3,crabA3);
grid on
box on
xlabel('H/km')
ylabel('ƫ����/^o')
legend('(c) \phi_c, �������SBR(\alpha_1=0^o, \eta_i=45^o)')
subplot(3,2,6)
plot(H/1e3,crabM3);
grid on
box on
xlabel('H/km')
ylabel('ƫ������')
% axis([0,10000,1,1.8])
legend('(d) \rho_c, �������SBR(\alpha_1=0^o, \eta_i=45^o)')