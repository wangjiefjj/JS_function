%% ͼ4.33����б�����ƫ���Ǻ�ƫ��������γ�ȵĹ�ϵ(H=506km)
clc;clear;close all
%% ���ƽ̨��������1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [506,1000,7000];
alpha1 = (0:80)/180*pi; %�������µ�����γ��
eta = (0:80)/180*pi; % fai���״��������н�
[X,Y] = meshgrid(eta,alpha1);
%% ���㺽ƫ��, ��ƫ����
crabA1 = real(fun_CrabAngle(X,Y,H(1))/pi*180);
crabM1 = (fun_CrabMagnitude(X,Y,H(1)));
%% figure
figure(1)
mesh(X/pi*180,Y/pi*180,crabA1);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\phi_c')
title('��б�����ƫ������γ�ȵĹ�ϵ(H=506km)')
%% figure
figure(2)
mesh(X/pi*180,Y/pi*180,crabM1);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\rho_c')
title('��б�����ƫ��������γ�ȵĹ�ϵ(H=506km)')

%% ���㺽ƫ��, ��ƫ����
crabA2 = real(fun_CrabAngle(X,Y,H(2))/pi*180);
crabM2 = (fun_CrabMagnitude(X,Y,H(2)));
%% figure
figure(3)
mesh(X/pi*180,Y/pi*180,crabA2);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\phi_c')
title('��б�����ƫ������γ�ȵĹ�ϵ(H=1000km)')
%% figure
figure(4)
mesh(X/pi*180,Y/pi*180,crabM2);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\rho_c')
title('��б�����ƫ��������γ�ȵĹ�ϵ(H=1000km)')


%% ���㺽ƫ��, ��ƫ����
crabA3 = real(fun_CrabAngle(X,Y,H(3))/pi*180);
crabM3 = (fun_CrabMagnitude(X,Y,H(3)));
%% figure
figure(5)
mesh(X/pi*180,Y/pi*180,crabA3);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\phi_c')
title('��б�����ƫ������γ�ȵĹ�ϵ(H=7000km)')
%% figure
figure(6)
mesh(X/pi*180,Y/pi*180,crabM3);
xlabel('\eta_i')
ylabel('\alpha1')
zlabel('\rho_c')
title('��б�����ƫ��������γ�ȵĹ�ϵ(H=7000km)')