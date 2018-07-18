%% ͼ4.38��ƫ���ǶԶ�����Ƶ��������Ӱ��,mesh ͼ����λ��-30:30
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                                  %��Ƶ Hz
c = 3e8;                                    %���� m/s
lambda = c/fo;                              %���� m
PRF = 500;                                  %�����ظ�Ƶ�� Hz
Tr = 1/PRF;                                 %�����ظ���� s
H = 700e3;                                    %���Ǹ߶� km
angle_Az = [0:1:180]./180*pi;              %��λ�� rad
alpha1 = 30/180*pi;                         %γ�� rad  
eta = 70/180*pi;                            %������ rad
R = 0:1e3:3000e3;                                 %�ؾ� km
for i = 1:length(angle_Az)
%% ��������ת������
    wd(i,:) = fun_Wd(H,R,angle_Az(i), Tr, lambda);
    %% ��������ת������
    wdr(i,:) = fun_Wd_Rotation(H,R,angle_Az(i), Tr, lambda, alpha1, eta);
end
%% figure
figure(1)
[X,Y] = meshgrid(R/1e3,angle_Az./pi.*180);
mesh(X,Y, (wdr))
hold on
mesh(X,Y, (wd))
% axis([0,3000,-15,60])
grid on
box on
xlabel('�ؾ�/km')
ylabel('��λ��/deg')
zlabel('������/Hz')
title(['ƫ���ǶԶ�����Ƶ��������Ӱ��(H=',num2str(H/1e3),'km)'])