%% R=500km的角度-多普勒域的理想杂波谱
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1.25e9;                                  %载频 Hz
C = 3e8;                                    %光速 m/s
lambda = C/fo;                              %波长 m
N = 32;                             %阵元个数
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 506;                            %SBR高度km
Vp = fun_Vp(H);                     %SBR速度m/s
d = 13.4;                           %归一化阵元间距
PRF = 500;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
c = 3e5;                            %光速km/m
beta = 19.47;
beta0 = beta * d;
% beta0 = 2*Vp*Tr/lambda/2;
alpha1 = 45/180*pi;                  %卫星纬度rad
eta = 90/180*pi;                    %卫星倾角rad
cosAz = linspace(-0.05,0.05,100);
RR = 500;                           %地距km 
EL = fun_ELAngle(H,RR)/180*pi;       %俯仰角
wd = -linspace(-1,1,200);                       %规1化多普勒
% Pt = zeros(length(wd),length(Az));
tic
for i_Az = 1:length(cosAz)                 %俯仰角
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %入射锥角
    %% 角度-多普勒域的理想杂波谱 
    %% 杂波协方差
    wdc = beta0*cmj;                         %杂波归1化多普勒
    b = exp(-1j*(0:M-1)*pi*wdc).';           %时间导向矢量
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %空间导向矢量 
    sc = kron(b,a);                          %杂波的导向矢量
    Rk = fun_SCMN(sc);                       %杂波协方差
%     Rk = Rk + eye(M*N);
%     iRk = inv(Rk);
    %% 在cos(Az(i_Az))方位角下的杂波功率普
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %时间导向矢量
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %空间导向矢量
        s = kron(b,a);                           %遍历的导向矢量
        Pt(i_wd,i_Az) = (abs(s'*Rk*s));          % 杂波功率普6.28    
    end
end
toc
%% figure
figure(1)
imagesc(cosAz,wd,abs((Pt)))
xlabel('cos(\theta)')
ylabel('归一化多普勒')
title('R=500km')
%% R=1400
RR = 1400;
EL = fun_ELAngle(H,RR)/180*pi;             %俯仰角
for i_Az = 1:length(cosAz)                 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %入射锥角
    %% 角度-多普勒域的理想杂波谱 
    %% 杂波协方差
    wdc = beta0*cmj;                         %杂波归1化多普勒
    b = exp(-1j*(0:M-1)*pi*wdc).';           %时间导向矢量
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %空间导向矢量 
    sc = kron(b,a);                          %杂波的导向矢量
    sc = awgn(sc,40);
    Rk = fun_SCMN(sc);                       %杂波协方差
    %% 在cos(Az(i_Az))方位角下的杂波功率普
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %时间导向矢量
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %空间导向矢量
        s = kron(b,a);                           %遍历的导向矢量
        Pt2(i_wd,i_Az) = (abs(s'*Rk*s));          % 杂波功率普6.28    
    end
end

%% figure
figure(2)
imagesc(cosAz,wd,abs((Pt2)))
xlabel('cos(\theta)')
ylabel('归一化多普勒')
title('R=1400km')