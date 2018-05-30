%% ͼ4.38��ƫ���ǶԶ�����Ƶ��������Ӱ��
clc;clear;close all
%% ���ƽ̨��������2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
fo = 40e9; %��ƵHz
c = 3e8;% m
lambda = c/fo;
PRF = 500;
Tr = 1/PRF;
H = 500;%km
angle_Az = [75,80,90]./180*pi; % ��λ��
cosAz = cos(angle_Az);
alpha1 = 30/180*pi;%γ��rad  
eta = 45/180*pi; %������rad
R = 0:3000; %�ؾ�km
%% ��������ת������
% wd1 = fun_Wd(H,R,angle_Az(1), Tr, lambda);
wd2 = fun_Wd(H,R,angle_Az(1), Tr, lambda);
wd3 = fun_Wd(H,R,angle_Az(2), Tr, lambda);
wd4 = fun_Wd(H,R,angle_Az(3), Tr, lambda);
%% ��������ת������
% wdr1 = fun_Wd_Rotation(H,R,angle_Az(1), Tr, lambda, alpha1, eta);
wdr2 = fun_Wd_Rotation(H,R,angle_Az(1), Tr, lambda, alpha1, eta);
wdr3 = fun_Wd_Rotation(H,R,angle_Az(2), Tr, lambda, alpha1, eta);
wdr4 = fun_Wd_Rotation(H,R,angle_Az(3), Tr, lambda, alpha1, eta);
%% figure
figure(1)
hold on
% plot(R,wd1,'r--','LineWidth',2);
% plot(R,wdr1,'r','LineWidth',2);

plot(R,wd2,'k--','LineWidth',2);
plot(R,wdr2,'k','LineWidth',2);

plot(R,wd3,'g--','LineWidth',2);
plot(R,wdr3,'g','LineWidth',2);

plot(R,wd4,'b--','LineWidth',2);
plot(R,wdr4,'b','LineWidth',2);

h_leg = legend('����ת\theta_{AZ}=75^{o}','������ת\theta_{AZ}=75^{o}',...
    '����ת\theta_{AZ}=80^{o}','������ת\theta_{AZ}=80^{o}',...
    '����ת\theta_{AZ}=90^{o}','������ת\theta_{AZ}=90^{o}');
set(h_leg,'Location','SouthEast')
axis([0,3000,-15,60])
grid on
box on
xlabel('�ؾ�R/km')
ylabel('������/Hz')
title('ͼ4.38 ƫ���ǶԶ�����Ƶ��������Ӱ��(H=500km)')
