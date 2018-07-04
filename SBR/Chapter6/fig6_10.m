%% 方向图
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                         %载频 Hz
C = 3e8;                            %光速 m/s
lambda = C/fo;                      %波长 m
N = 16;                             %阵元个数
P = 4;
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 1250;                           %SBR高度km
Vp = fun_Vp(H);                     %SBR速度m/s
d = 1;                            %归一化阵元间距
% d = 8.3;                          %归一化阵元间距
gamma = 5;                         %发射接收空间比率
PRF = 400;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
c = 3e5;                            %光速km/m
% beta = 19.47;
beta = Vp*4/PRF/lambda;
beta0 = beta * d;
%% 方向图 
theta0_AZ = 10/180*pi;                 %方位向指向
% theta0_EL = 10/180*pi;                 %俯仰向指向
theta_AZ_du = (-100:1:100);
theta_AZ=(theta_AZ_du)./180*pi;           %方位向搜索范围和步长
theta_EL=(0:40)./180*pi;                %俯仰向搜索范围和步长
% w = ones(N*P,1);
% window = taylorwin(N*P);
% s = s.*window;
for i_EL = 1:length(theta_EL)
    i_EL
    c_s = cos(theta_EL(i_EL))*sin(theta0_AZ);   %信号锥角 sin(theta0_EL)*
    a = exp(-1j*(0:N-1)*pi*d*c_s).';   %
    b = exp(-1j*(0:P-1)*pi*gamma*d*c_s).';
    w = kron(b,a);                  %权值
    w = w/sqrt(w'*w);               %归一化
    cmj = cos(theta_EL(i_EL))*sin(theta_AZ); 
    for i_AZ = 1:length(theta_AZ)
        a = exp(-1j*(0:N-1)*pi*d*cmj(i_AZ)).';               %接收空间导向矢量 
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj(i_AZ)).';         %发射空间导向矢量
        s = kron(b,a);                                       %信号
        s = s/sqrt(s'*s);                                 %归一化
%         w = w.*window;
        A(i_AZ,i_EL) = w'*s;                          %方向图(每一行是一个方位角)       
    end
end
figure()
plot(theta_AZ_du,10*log10(abs(A(:,1))))
[hang,lie] = size(A);
if lie>1
    figure()
    mesh(10*log10(abs(A)))
end
% axis([-100,100,-1e-7,0])