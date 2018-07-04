%% 角度-多普勒域的理想杂波谱
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                                  %载频 Hz
C = 3e8;                                    %光速 m/s
lambda = C/fo;                              %波长 m
N = 16;                             %阵元个数
P = 4;
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 7000;                            %SBR高度km
Vp = fun_Vp(H);                     %SBR速度m/s
d = 63.9;                           %归一化阵元间距
% d = 8.3;                           %归一化阵元间距
gamma = 5;                              %发射接收空间比率
PRF = 400;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
c = 3e5;                            %光速km/m
% beta = 19.47;
beta = Vp*4/PRF/lambda;
beta0 = beta * d;

% beta0 = 2*Vp*Tr/lambda/2;
alpha1 = 45/180*pi;                  %卫星纬度rad
eta = 90/180*pi;                    %卫星倾角rad
Az = linspace(-27/180*pi,27/180*pi,100);
% Az = linspace(0.4,0.42,100);
cosAz = cos(Az);
% cosAz = linspace(-0.05,0.05,100);
%% R=500
RR = 500;                           %地距km 
EL = fun_ELAngle(H,RR)/180*pi;       %俯仰角
wd = -linspace(-1,1,200);                       %规1化多普勒
% Pt = zeros(length(wd),length(Az));
tic
for i_Az = 1:length(cosAz)                 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %入射锥角
    %% 角度-多普勒域的理想杂波谱 
    %% 杂波协方差
    wdc = beta0*cmj;                         %杂波归1化多普勒
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %接收空间导向矢量 
    b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).';         %发射空间导向矢量  
    c = exp(-1j*(0:M-1)*pi*wdc).';           %时间导向矢量 
    sc = kron(c,kron(b,a));             %杂波的导向矢量
    sc = awgn(sc,20);
    Rk = fun_SCMN(sc);                       %杂波协方差
    %% 在cos(Az(i_Az))方位角下的杂波功率普
    for i_wd = 1:length(wd)
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %接收空间导向矢量
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).'./sqrt(P);         %发射空间导向矢量
        c = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %时间导向矢量
        s = kron(c,kron(b,a));                           %遍历的导向矢量
        Pt(i_wd,i_Az) = (abs(s'*Rk*s));          % 杂波功率普6.28    
    end
end
toc
%% figure
figure(1)
imagesc(Az/pi*180,wd,abs((Pt)))
xlabel('\theta_{AZ}/deg')
ylabel('归一化多普勒')
title(['H=',num2str(H),'km, R=',num2str(RR),'km'])
%% R=1400
RR = 1400;
EL = fun_ELAngle(H,RR)/180*pi;             %俯仰角
for i_Az = 1:length(cosAz)                 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %入射锥角
    %% 角度-多普勒域的理想杂波谱 
    %% 杂波协方差
    wdc = beta0*cmj;                         %杂波归1化多普勒
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %接收空间导向矢量 
    b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).';         %发射空间导向矢量  
    c = exp(-1j*(0:M-1)*pi*wdc).';           %时间导向矢量 
    sc = kron(c,kron(b,a));             %杂波的导向矢量
    sc = awgn(sc,20);
    Rk = fun_SCMN(sc);                       %杂波协方差
    %% 在cos(Az(i_Az))方位角下的杂波功率普
    for i_wd = 1:length(wd)
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %接收空间导向矢量
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).'./sqrt(P);         %发射空间导向矢量
        c = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %时间导向矢量
        s = kron(c,kron(b,a));                           %遍历的导向矢量
        Pt2(i_wd,i_Az) = (abs(s'*Rk*s));          % 杂波功率普6.28    
    end
end

%% figure
figure(2)
imagesc(Az/pi*180,wd,abs((Pt2)))
xlabel('\theta_{AZ}/deg')
ylabel('归一化多普勒')
title(['H=',num2str(H),'km, R=',num2str(RR),'km'])