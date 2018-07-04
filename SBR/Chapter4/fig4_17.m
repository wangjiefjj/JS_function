%% ͼ4.17��PFR���ؾ��ģ������Ĺ�ϵ
clc;clear;close all
%% ��������1
PRF1 = 400; %Hz
Tr1 = 1/PRF1; %s
H = [700e3, 1000e3, 2000e3, 7000e3]; %m
c = 3e8; %m
%% �����ӵؽ�
[Graze1,R1,~] = fun_GrazeAngle(H(1));
[Graze2,R2,~] = fun_GrazeAngle(H(2));
[Graze3,R3,~] = fun_GrazeAngle(H(3));
[Graze4,R4,~] = fun_GrazeAngle(H(4));
Graze1=Graze1/180*pi;
Graze2=Graze2/180*pi;
Graze3=Graze3/180*pi;
Graze4=Graze4/180*pi;
%% �������ģ������ ʽ 4.35
dR11 =  c*Tr1/2 *sec(Graze1);
dR21 =  c*Tr1/2 *sec(Graze2);
dR31 =  c*Tr1/2 *sec(Graze3);
dR41 =  c*Tr1/2 *sec(Graze4);
%% figure
figure(1)
% suptitle('�ܱ���')
set(gcf,'Position',[700 0 600 1000])
subplot(2,1,1)
hold on
plot(R1(500:end),dR11(500:end),'r','LineWidth',2);
plot(R2(500:end),dR21(500:end),'k','LineWidth',2);
plot(R3(500:end),dR31(500:end),'g','LineWidth',2);
plot(R4(500:end),dR41(500:end),'b','LineWidth',2);
h_leg = legend('H=700km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
axis([500,4000,300,1000])
grid on
box on
xlabel('�ؾ�R/km')
ylabel('\Delta_R/km')
title('ͼ4.17��a�����ģ������\Delta_R��ؾ�Ĺ�ϵ, PRF=500Hz')

%% ��������2
PRF2 = 2000; %Hz
Tr2 = 1/PRF2; %s
%% �������ģ������ ʽ 4.35
dR12 =  c*Tr2/2 *sec(Graze1);
dR22 =  c*Tr2/2 *sec(Graze2);
dR32 =  c*Tr2/2 *sec(Graze3);
dR42 =  c*Tr2/2 *sec(Graze4);
%% figure
subplot(2,1,2)
hold on
plot(R1(500:end),dR12(500:end),'r','LineWidth',2);
plot(R2(500:end),dR22(500:end),'k','LineWidth',2);
plot(R3(500:end),dR32(500:end),'g','LineWidth',2);
plot(R4(500:end),dR42(500:end),'b','LineWidth',2);
h_leg = legend('H=700km','H=1000km','H=2000km','H=7000km');
set(h_leg,'Location','SouthEast')
axis([500,4000,50,600])
grid on
box on
xlabel('�ؾ�R/km')
ylabel('\Delta_R/km')
title('ͼ4.17��b�����ģ������\Delta_R��ؾ�Ĺ�ϵ, PRF=2000Hz')
