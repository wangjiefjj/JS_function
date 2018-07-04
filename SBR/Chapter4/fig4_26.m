%% ͼ4.25�����뼰�Ӳ�б��(\beta_0=6)����Ķ�������չ
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1.25e9; %��ƵHz
c = 3e8;% m
lambda = c/fo;
H = 500e3;%m
Vp = fun_Vp(H);
angle_Az = (0:180)./180*pi; % ��λ��
R = [1000e3,500e3]; %�ؾ�km
beta0 = 6;
Tr =  lambda/2*beta0/2/Vp;
wd1 = mod(fun_Wd_beta0(H,R(1),angle_Az,beta0),2)-1;
wd2 = mod(fun_Wd_beta0(H,R(2),angle_Az,beta0),2)-1;
%% figure
figure(1)
hold on
plot(cos(angle_Az),wd1,'r','LineWidth',2);
plot(cos(angle_Az),wd2,'k','LineWidth',2);
h_leg = legend('R=1000km','R=500km');
set(h_leg,'Location','SouthEast')
% axis([100,2000,0,50])
grid on
box on
xlabel('Azimuth')
ylabel('������/Hz')
title('ͼ4.26 ���뼰�Ӳ�б��(\beta_0=6)����Ķ�������չ')