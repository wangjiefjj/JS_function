%% ͼ4.38��ƫ���ǶԶ�����Ƶ��������Ӱ��
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                                  %��Ƶ Hz
c = 3e8;                                    %���� m/s
lambda = c/fo;                              %���� m
PRF = 500;                                  %�����ظ�Ƶ�� Hz
Tr = 1/PRF;                                 %�����ظ���� s
H = 700e3;                                    %���Ǹ߶� km
angle_Az = [85,90,95]./180*pi;              %��λ�� rad
alpha1 = 30/180*pi;                         %γ�� rad  
eta = 70/180*pi;                            %������ rad
R = 0:1e3:3000e3;                           %�ؾ� m
%% ��������ת������
wd1 = fun_Wd(H,R,angle_Az(1), Tr, lambda);
wd2 = fun_Wd(H,R,angle_Az(2), Tr, lambda);
wd3 = fun_Wd(H,R,angle_Az(3), Tr, lambda);
%% ��������ת������
wdr1 = fun_Wd_Rotation(H,R,angle_Az(1), Tr, lambda, alpha1, eta);
wdr3 = fun_Wd_Rotation(H,R,angle_Az(2), Tr, lambda, alpha1, eta);
wdr4 = fun_Wd_Rotation(H,R,angle_Az(3), Tr, lambda, alpha1, eta);
%% figure
figure(1)
hold on
R = R/1e3;
plot(R,wd1,'r--','LineWidth',2);
plot(R,wdr1,'r','LineWidth',2);

plot(R,wd2,'g--','LineWidth',2);
plot(R,wdr3,'g','LineWidth',2);

plot(R,wd3,'b--','LineWidth',2);
plot(R,wdr4,'b','LineWidth',2);

h_leg = legend('����ת\theta_{AZ}=85^{o}','������ת\theta_{AZ}=85^{o}',...
    '����ת\theta_{AZ}=90^{o}','������ת\theta_{AZ}=90^{o}',...
    '����ת\theta_{AZ}=95^{o}','������ת\theta_{AZ}=95^{o}');
set(h_leg,'Location','SouthEast')
% axis([0,3000,-15,60])
grid on
box on
xlabel('�ؾ�R/km')
ylabel('������/Hz')
title(['ƫ���ǶԶ�����Ƶ��������Ӱ��(H=',num2str(H/1e3),'km)'])
