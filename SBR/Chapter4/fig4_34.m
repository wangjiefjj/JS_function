%% ͼ4.34����ͬ�������£���Ƭ����ƫ��������γ�ȵĹ�ϵ��H=506��
clc;clear;close all
%% ���ƽ̨��������1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [506,1000,7000];
eta = [90;60;30;-1]; %������rad
L = 200;
%% ����
Alpha1=[];
for i = 1:length(eta)
    alpha1 = linspace(-eta(i),eta(i),L); %γ��deg
    Alpha1 = [Alpha1;alpha1];
    %% ƫ����
    CrabA1(i,:) = real(fun_CrabAngle(alpha1/180*pi, eta(i)/180*pi, H(1)))/pi*180;
    %% ƫ������
    CrabM1(i,:) = real(fun_CrabMagnitude(alpha1/180*pi, eta(i)/180*pi, H(1)));

end
CrabA1 = [CrabA1;-CrabA1];
Alpha1t = [Alpha1;Alpha1];
Eta = repmat(eta,[1,L]);
Etat = repmat(eta,[2,L]);
%% figure
figure(1)
contour(Alpha1t,CrabA1,Etat,[90,60,30,0],'ShowText','on','LineWidth',5)
grid on 
xlabel('γ��\alpha_1/^o')
ylabel('ƫ����\phi_c/^o')
axis([-100,100,-4,4])
%% figure
figure(2)
contour(Alpha1,CrabM1,Eta,[90,60,30,0],'ShowText','on','LineWidth',5)
grid on 
xlabel('γ��\alpha_1/^o')
ylabel('ƫ������\rho_c/^o')
axis([-100,100,0.93,1.01])
