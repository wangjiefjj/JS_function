%% 图6.6，无地球自转的 方位=89.5度使得2D波束形成时的杂波谱
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;                             %阵元个数
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 506;                            %SBR高度km
Vp = fun_Vp(H);                     %SBR速度
d = 13.4;                           %归一化阵元间距
PRF = 500;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
c = 3e5;                            %光速km/m
beta = 19.47;
beta0 = beta * d;
% alpha1 = 30/180*pi;                  %卫星纬度
% eta = 45/180*pi;                    %卫星倾角
%% 目标区域
RR = 400:10:1100;                    %地距km 
Az = 89.5/180*pi;                    %方位角 
EL = fun_ELAngle(H,RR)/180*pi;       %俯仰角
c = sin(EL)*cos(Az);                 %入射锥角
L = length(c);                       %距离单元

% c = -1:0.01:1;
% Az = 89.5/180*pi;                      %方位角 
% L = length(c);                       %距离单元
% RR = linspace(400,1100,L);   

%% 目标导向矢量计算
wdc = beta0*c/PRF;                 %规1化多普勒
S = zeros(M*N,L);
% w = hamming(L);
%杂波生成
v = 1;
for i = 1:length(c)
    b = exp(-1j*(0:M-1)*pi*wdc(i)).';         %时间导向矢量
    a = exp(-1j*(0:N-1)*pi*d*c(i)).';         %空间导向矢量 
    S(:,i) = gamrnd(v,1/v,1,1)*kron(b,a);     %杂波的导向矢量
%     X(:,i) = awgn(S(:,i),0);
end
R = fun_SCMN(S);
invR = inv(R);
wd = -1:0.01:1;
for i_c = 1:length(c)
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).';         %时间导向矢量
        a = exp(-1j*(0:N-1)*pi*d*c(i_c)).';        %空间导向矢量 
        s = kron(b,a);
        P(i_c,i_wd)= s'*R*s;
    end
end
%% figure
figure(1)
[X,Y]=meshgrid(wd,RR);
mesh(X,Y,abs(fftshift(P)))
xlabel('归一化多普勒')
ylabel('距离/km')
