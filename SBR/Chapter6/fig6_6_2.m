%% 图6.6， 方位=89.5度使得2D波束形成时的杂波谱,可以存在自转和不存在自转都能画
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9600e6; %450e6   9600e6                 %载频 Hz
C = 299792458;                            %光速 m/s
lambda = C/fo;                      %波长 m
Nr = 16;                             %接收阵元个数
Nt = 4;                              %发射阵元个数
Np = 16;                             %相干脉冲数
D = Nt*Nr*Np;                          %自由度 
Re = 6373e3;                          %地球半径m
H = 700e3;  %700e3  %9e3                 %SBR高度 m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR速度m/s
end                
d = 60;%                           %归一化阵元间距
gamma = 16;                         %发射接收空间比率
PRF = 400;                         %脉冲重复频率Hz，500~2000
Tr =1/PRF;                         %脉冲重复间隔s
% beta = 19.47;
beta = Vp*2/PRF/(d*lambda/2);      %%混叠系数 
alpha1 = 20/180*pi;                  %卫星纬度
eta = 90/180*pi;                    %卫星倾角

%% 偏航角幅度
CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度

%% 目标区域
RR = 1000e3:10e3:3400e3;                   %地距m 
Rs = fun_R2Rs(H,RR);
Graze = fun_GrazeAngle(H, R, Rs)/180*pi;
dR =  c*Tr/2 *sec(Graze);   
EL_du = fun_ELAngle(H,RR);           %波束俯仰角
EL = EL_du/180*pi;
Num = 360;                          %方位分块
CNR = 40;                   %杂噪比dB
wd = linspace(-0.5,0.5,1000);       %虚拟规1化多普勒
fsp = 0;
Pt = zeros(length(RR),length(wd));
Ptr = zeros(length(RR),length(wd));
Rk = 
for i = 1:length(EL)
    i
    Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma);%,CrabA,CrabM
    Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma,CrabA,CrabM); 
    for i_wd = 1:length(wd)
        a = exp(1j*(0:Nr-1)*2*pi*fsp).';               %接收空间导向矢量
        b = exp(1j*(0:Nt-1)*2*pi*gamma*fsp).';         %发射空间导向矢量
        c = exp(1j*(0:Np-1)*2*pi*wd(i_wd)).';            %时间导向矢量
        s = kron(c,kron(b,a));                                  %遍历的导向矢量
        Pt(i,i_wd)= (abs(s'*Rk*s));                  % 匹配滤波SINR 6.27
        Ptr(i,i_wd)= (abs(s'*Rk_r*s));                  % 匹配滤波SINR 6.27
    end
end

%% figure
figure(1)
[X,Y]=meshgrid(wd,RR);
mesh(X,Y,(abs((Pt))))%10*log10
xlabel('归一化多普勒')
ylabel('距离/km')
view(-0,90)
title(['R-D图，H=',num2str(H),'m'])
figure(2)
[X,Y]=meshgrid(wd,RR);
imagesc(wd,RR/1e3,(abs((Pt))))
xlabel('归一化多普勒')
ylabel('距离/km')
view(-0,-90)
title(['R-D图，H=',num2str(H),'m'])
figure(3)
imagesc(wd,RR/1e3,(abs((Ptr))))
xlabel('归一化多普勒')
ylabel('距离/km')
view(-0,-90)
title(['R-D图带自转，H=',num2str(H),'m'])