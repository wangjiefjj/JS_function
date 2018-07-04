%% 接收波束形成
%% 方向图
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                         %载频 Hz
C = 3e8;                            %光速 m/s
lambda = C/fo;                      %波长 m
Nrx = 16;                           %接收方位向阵元个数
Nry = 2;                            %接收俯仰向阵元个数
P = 4;                              %发射通道
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 1250;                           %SBR高度km
Vp = fun_Vp(H);                     %SBR速度m/s
drx = 63.9;                         %接收x轴归一化阵元间距
dry = 16;                           %接收y轴归一化阵元间距
% d = 8.3;                          %归一化阵元间距
gamma = 16;                          %发射接收空间比率(x方向)
PRF = 400;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
c = 3e5;                            %光速km/m
% beta = 19.47;
% beta = Vp*4/PRF/lambda;
% beta0 = beta * d;
%% 期望信号
theta0_AZ = 50/180*pi;                 %方位向指向
theta0_EL = 20/180*pi;                 %俯仰向指向

theta_AZ_du = (-100:1:100);
theta_EL_du = (-100:1:100);
theta_AZ=(theta_AZ_du)./180*pi;         %方位向搜索范围和步长
theta_EL=(theta_EL_du)./180*pi;                %俯仰向搜索范围和步长
%% 权值
c_sx = sin(theta0_AZ);   %信号锥角 sin(theta0_EL)*
c_sy = sin(theta0_AZ)*sin(theta0_EL);   %信号锥角 sin(theta0_EL)*
ax = exp(-1j*(0:Nrx-1)*pi*drx*c_sx).';   %
ay = exp(-1j*(0:Nry-1)*pi*dry*c_sy).';   %
a = kron(ax,ay);
w = a;
% b = exp(-1j*(0:P-1)*pi*gamma*drx*c_sx).';
% w = kron(b,a);                  %权值
% w = w/sqrt(w'*w);               %归一化
%% 方向图
for i_EL = 1:length(theta_EL)
    i_EL
    cx =  sin(theta_AZ);
    cy =  sin(theta_AZ)*sin(theta_EL(i_EL)); %sin(theta_EL(i_EL))* sin(theta_AZ)
    for i_AZ = 1:length(theta_AZ)        
        ax = exp(-1j*(0:Nrx-1)*pi*drx*cx(i_AZ)).';               %接收空间导向矢量x轴 
        ay = exp(-1j*(0:Nry-1)*pi*dry*cy(i_AZ)).';   %
        a = kron(ax,ay);
        s = a;
%         b = exp(-1j*(0:P-1)*pi*gamma*drx*cx(i_AZ)).';         %发射空间导向矢量
%         s = kron(b,a);                                       %信号
%         s = s/sqrt(s'*s);                                 %归一化
%         w = w.*window;
        A(i_AZ,i_EL) = w'*s;                          %方向图(每一行是一个方位角)       
    end
end
figure()
plot(theta_AZ_du,10*log10(abs(A(:,1))))
figure()
plot(theta_EL_du,10*log10(abs(A(1,:))))
[hang,lie] = size(A);
if lie>1
    figure()
    mesh(10*log10(abs(A)))
end
% axis([-100,100,-1e-7,0])