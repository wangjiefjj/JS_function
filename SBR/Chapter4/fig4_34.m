%% 图4.34，不同轨道倾角下，航片角与偏航幅度与纬度的关系（H=506）
clc;clear;close all
%% 天基平台参数设置1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [506,1000,7000];
eta = [90;60;30;-1]; %轨道倾角rad
L = 200;
%% 计算
Alpha1=[];
for i = 1:length(eta)
    alpha1 = linspace(-eta(i),eta(i),L); %纬度deg
    Alpha1 = [Alpha1;alpha1];
    %% 偏航角
    CrabA1(i,:) = real(fun_CrabAngle(alpha1/180*pi, eta(i)/180*pi, H(1)))/pi*180;
    %% 偏航幅度
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
xlabel('纬度\alpha_1/^o')
ylabel('偏航角\phi_c/^o')
axis([-100,100,-4,4])
%% figure
figure(2)
contour(Alpha1,CrabM1,Eta,[90,60,30,0],'ShowText','on','LineWidth',5)
grid on 
xlabel('纬度\alpha_1/^o')
ylabel('偏航幅度\rho_c/^o')
axis([-100,100,0.93,1.01])
