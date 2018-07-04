%% 图6.6， 方位=89.5度使得2D波束形成时的杂波谱,可以存在自转和不存在自转都能画
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                                  %载频 Hz
C = 3e8;                                    %光速 m/s
lambda = C/fo;                              %波长 m
N = 16;                             %阵元个数
P = 4;
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 700e3;                            %SBR高度km
Vp = fun_Vp(H);                     %SBR速度m/s
d = 1;                           %归一化阵元间距
% d = 8.3;                           %归一化阵元间距
gamma = 5;                              %发射接收空间比率
PRF = 400;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
c = 3e5;                            %光速km/m
% beta = 19.47;
beta = Vp*2/PRF/(d*2*lambda);
beta0 = beta * d;

alpha1 = 45/180*pi;                  %卫星纬度
eta = 90/180*pi;                    %卫星倾角

%% 目标区域
RR = 800e3:10e3:2400e3;                    %地距km 
wd = -1:0.01:1;                      %仿真的多普勒范围
PP = zeros(length(RR),length(wd));
for Az = 0:0.1:0                   %方位角/deg
    Az
Az = Az/180*pi;
EL = fun_ELAngle(H,RR)/180*pi;       %俯仰角
CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度
cmj = sin(EL)*cos(Az);                     %入射锥角
% cmj = CrabM*sin(EL)*cos(Az+CrabA);       %入射锥角(考虑自转)
L = length(cmj);                       %距离单元

%% 杂波RD图
wdc = beta0*cmj;                 %杂波归1化多普勒
S = zeros(P*M*N,L);
tic
Pt = zeros(length(cmj),length(wd));
for i = 1:length(cmj)
    i
    a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %接收空间导向矢量 
    b = exp(-1j*(0:P-1)*pi*gamma*d*cmj(i)).';         %发射空间导向矢量  
    c = exp(-1j*(0:M-1)*pi*wdc(i)).';           %时间导向矢量 
    S(:,i) = kron(c,kron(b,a));                         %杂波的导向矢量
    X(:,i) = awgn(S(:,i),20);
    R = fun_SCMN(X(:,i));                       %不同距离杂波协方差
    for i_wd = 1:length(wd)
        a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %接收空间导向矢量 
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj(i)).';         %发射空间导向矢量  
        c = exp(-1j*(0:M-1)*pi*wd(i_wd)).';         %时间导向矢量    
        s = kron(c,kron(b,a));                              %遍历的导向矢量
        Pt(i,i_wd)= (abs(s'*R*s));                  % 匹配滤波SINR 6.27
    end
end
toc
PP = PP+Pt;
end
%% figure
% figure(1)
% [X,Y]=meshgrid(wd,RR);
% mesh(X,Y,abs((P)))
% xlabel('归一化多普勒')
% ylabel('距离/km')
% view(-0,90)
% title(['R-D图，H=',num2str(H),'km'])
figure(2)
[X,Y]=meshgrid(wd,RR);
imagesc(wd,RR,abs((PP)))
xlabel('归一化多普勒')
ylabel('距离/km')
view(-0,-90)
title(['R-D图，H=',num2str(H),'km'])