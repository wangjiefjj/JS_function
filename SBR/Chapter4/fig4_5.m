%% ͼ4.5 ���ؾ࣬���������߶ȵĹ�ϵ
clc;clear;close all
%% ���ƽ̨����
H = 100:10000; %% ����״�����µ�ľ���,km
%% �������ؾ�
Rmax = fun_Rmax(H);
%% figure
figure(1)
plot(H,Rmax ,'LineWidth',2)
grid on
box on
xlabel('�߶�H/km')
ylabel('���ؾ�Rmax/km')
title('ͼ4.5 ���ؾ���߶ȵĹ�ϵ')
%% �����������
ELmax = fun_ELmax(H);
%% figure
figure(2)
plot(H,ELmax ,'LineWidth',2)
grid on
box on
xlabel('�߶�H/km')
ylabel('�������\theta_{EL,max}/km')
title('ͼ4.5 ���������߶ȵĹ�ϵ')
