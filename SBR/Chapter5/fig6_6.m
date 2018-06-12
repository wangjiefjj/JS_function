%% 图6.6， 方位=89.5度使得2D波束形成时的杂波谱,可以存在自转和不存在自转都能画
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1.25e9;                                  %载频 Hz
C = 3e8;                                    %光速 m/s
lambda = C/fo;                              %波长 m
N = 12;                             %阵元个数
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
alpha1 = 45/180*pi;                  %卫星纬度
eta = 90/180*pi;                    %卫星倾角

%% 目标区域
RR = 400:10:1100;                    %地距km 
wd = -1:0.001:1;                      %仿真的多普勒范围
P = zeros(length(RR),length(wd));
for Az = 89.5:0.1:89.5                   %方位角 
    Az
Az = Az/180*pi;
EL = fun_ELAngle(H,RR)/180*pi;       %俯仰角
CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度
% cmj = sin(EL)*cos(Az);                 %入射锥角
cmj = CrabM*sin(EL)*cos(Az+CrabA);       %入射锥角(考虑自转)
L = length(cmj);                       %距离单元

%% 杂波RD图
wdc = beta0*cmj;                 %杂波归1化多普勒
S = zeros(M*N,L);
tic
Pt = zeros(length(cmj),length(wd));
for i = 1:length(cmj)
    b = exp(-1j*(0:M-1)*pi*wdc(i)).';           %时间导向矢量
    a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %空间导向矢量 
    S(:,i) = kron(b,a);                         %杂波的导向矢量
    X(:,i) = awgn(S(:,i),0);
    R = fun_SCMN(X(:,i));                       %不同距离杂波协方差
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).';         %时间导向矢量
        a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %空间导向矢量 
        s = kron(b,a);                              %遍历的导向矢量
        Pt(i,i_wd)= (abs(s'*R*s));                  % 匹配滤波SINR 6.27
    end
end
toc
P = P+Pt;
end
%% figure
figure(1)
[X,Y]=meshgrid(wd,RR);
mesh(X,Y,abs((P)))
xlabel('归一化多普勒')
ylabel('距离/km')
view(-0,90)