%% R=500km的角度-多普勒域的理想杂波谱
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.26e9;                                  %载频 Hz
C = 3e8;                                    %光速 m/s
lambda = C/fo;                              %波长 m
N = 32;                             %阵元个数
M = 16;                             %相干脉冲数
Re = 6373;                          %地球半径km
H = 1250;                            %SBR高度km
Vp = fun_Vp(H);                     %SBR速度m/s
d = 63.9;                           %归一化阵元间距
% d = 10;                            %归一化阵元间距
PRF = 500;                          %脉冲重复频率Hz，500~2000
Tr =1/PRF;                          %脉冲重复间隔s
beta =  Vp*2/PRF/lambda;             %
beta0 = beta * d;

alpha1 = 45/180*pi;                  %卫星纬度rad
eta = 90/180*pi;                     %卫星倾角rad

%% 杂波协方差
RR = 500;                           %地距km 
EL = fun_ELAngle(H,RR)/180*pi;       %俯仰角
Rk = fun_GenerateR(1, N, M, 256, 40, EL, beta0, d);
% Az = linspace(0,180,256)/180*pi;     %方位角deg
% % cosAz = cos(Az);
% %% 锥角和杂波多普勒
% CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
% CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度
% % cmj = sin(EL).*cos(AZ);         %入射锥角
% cosAz = cos(Az+CrabA);



% tic
% %% 杂波协方差
% Rk = zeros(M*N,M*N);
% CNR = 40;                                  %杂噪比
% for i_Az = 1:length(cosAz)                 %俯仰角
%     i_Az
%     cmj = sin(EL)*CrabM*cosAz(i_Az);               %入射锥角
%     wdc = beta0*cmj;                         %杂波归1化多普勒
%     b = exp(-1j*(0:M-1)*pi*wdc).';           %时间导向矢量
%     a = exp(-1j*(0:N-1)*pi*d*cmj).';         %空间导向矢量 
%     sc = kron(b,a);                          %杂波的导向矢量
%     Rk = Rk + fun_SCMN(sc);                  %杂波协方差
% end
% Rk = Rk/max(max(abs(Rk)));
% Rk = Rk*10^(CNR/10) + eye(size(Rk));
% iRk = inv(Rk);

%% 角度-多普勒域的理想杂波谱
% Az = linspace(0,180,256)/180*pi;     %方位角deg
% cosAz = cos(Az);
cosAz = linspace(-1,1,256);
wd = -linspace(-1,1,200);                       %规1化多普勒
for i_Az = 1:length(cosAz) 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);               %入射锥角
    for i_wd = 1:length(wd)
         b = exp(-1j*(0:M-1)*pi*wd(i_wd)).';           %时间导向矢量
         a = exp(-1j*(0:N-1)*pi*d*cmj).';         %空间导向矢量
         s = kron(b,a);                           %遍历的导向矢量
         Pt(i_wd,i_Az) =(abs(s'*Rk*s));          % 杂波功率普6.28    
    end
end
%% figure
figure(1)
[X,Y]=meshgrid(cosAz,wd);
mesh(X,Y,((Pt)))
xlabel('cos(\theta)')
ylabel('归一化多普勒')
view(-0,90)

figure(2)
imagesc(cosAz,wd,((Pt)))
xlabel('cos(\theta)')
ylabel('归一化多普勒')