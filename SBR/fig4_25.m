%% ͼ4.25�����뷽λ��ȶ�������
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 4.5e9; %��ƵHz
c = 3e8;% m
lambda = c/fo;
PRF = 500;
Tr = 1/PRF;
H = 500;%km
angle_Az = (0:180)./180*pi; % ��λ��
cosAz = cos(angle_Az);
R = 0:3000; %�ؾ�km
%% ���������
for i = 1:length(angle_Az)
    wd(i,:) = fun_Wd(H,R,angle_Az(i), Tr, lambda);
end
%% figure
figure(1)
contour(R,angle_Az/pi*180,wd,20)
xlabel('�ؾ�R/km')
ylabel('��λ��/^{o}')
title('ͼ4.25 ���뷽λ��ȶ�������')
%% figure
figure(2)
contour(R,cosAz,wd,20)
xlabel('�ؾ�R/km')
ylabel('cos(\theta_{AZ})')
title('ͼ4.25 ���뷽λ��ȶ�������')
%% figure
figure(3)
[X,Y] = meshgrid(R,angle_Az/pi*180);
mesh(X,Y,wd)
xlabel('�ؾ�R/km')
ylabel('��λ��/^{o}')
zlabel('������/Hz')
title('ͼ4.25 ���뷽λ��ȶ�������')
