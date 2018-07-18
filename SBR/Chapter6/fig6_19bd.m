%% 图6.19， 方位角90的RD图，有距离模糊有自转
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9600e6; %450e6   9600e6                 %载频 Hz
c = 299792458;                            %光速 m/s
lambda = c/fo;                      %波长 m
Nr = 32;                             %接收阵元个数
Nt = 1;                              %发射阵元个数
Np = 32;                             %相干脉冲数
D = Nt*Nr*Np;                          %自由度 
Re = 6373e3;                          %地球半径m
H = 700e3;  %700e3  %9e3                 %SBR高度 m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR速度m/s
end                
d = 1;%                           %归一化阵元间距
gamma = 8;                         %发射接收空间比率
fr = 500;                         %脉冲重复频率Hz，500~2000
Tr =1/fr;                         %脉冲重复间隔s
% beta = 19.47;
beta = Vp*2/fr/(d*lambda/2);      %%混叠系数 
alpha1 = 20/180*pi;                  %卫星纬度
eta = 90/180*pi;                    %卫星倾角

%% 偏航角幅度
CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度

%% 目标区域
RR = 1000e3:10e3:2400e3;                   %地距m 
Rs = fun_R2Rs(H,RR);
% Graze = fun_GrazeAngle(H, RR, Rs)/180*pi;
% dR =  C*Tr/2 *sec(Graze);
% %% 模糊距离和俯仰角
% EL_du(1,:) = fun_ELAngle(H,RR);           %波束俯仰角
% EL(1,:) = EL_du(1,:)/180*pi;
Rsmax = fun_Rsmax(H);                       %最大探测斜距m
Rsu=c/(2*fr);                               %最大无模糊斜距m
% for i = 2:2:4
%     RRt = RR + (i-1)*dR;
%     EL_du(i,:) = fun_ELAngle(H,RRt);           %波束俯仰角
%     EL(i,:) = EL_du(i,:)/180*pi;
%     RRt = RR - (i-1)*dR;
%     EL_du(i+1,:) = fun_ELAngle(H,RRt);           %波束俯仰角
%     EL(i+1,:) = EL_du(i,:)/180*pi;   
% end
%% 杂波
[ ~,El0,d,lambda,fr,Np,Nr ] = fun_JWR( 0,0);
Num = 180;                          %方位分块
CNR = 40;                   %杂噪比dB
wd = linspace(-0.5,0.5,100);       %虚拟规1化多普勒
cmj = wd./beta/(d*lambda^2);
fsp = 0;%;d*cmj;
Pt = zeros(length(RR),length(wd));
Ptr = zeros(length(RR),length(wd));
for i = 1:length(RR)
    i
%     Nk1 = floor((Rsmax-Rs(i))/Rsu);                %前向模糊数
%     Nk2 = floor((Rs(i)-H)/Rsu);                    %后向模糊数
%     Rsk1 = Rs(i) + (0:Nk1)*Rsu;                    %前向模糊斜距m
%     Rsk2 = Rs(i) - (Nk2:-1:1)*Rsu;                 %后向模糊斜距m
%     Rsk = [Rsk2,Rsk1];%Rs(i);%                 %模糊斜距m
%     RRk = fun_Rs2R(H,Rsk);                       %模糊地距m
%     elk = fun_ELAngle(H,RRk)/180*pi;             %各模糊距离俯仰角rad
%     Rk = zeros(D,D);
%     Rk_r = zeros(D,D);
%     for j = 1:length(elk)
%         Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(j), beta, d,gamma, lambda);%,CrabA,CrabM
%         Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(j), beta, d,gamma,lambda,CrabA,CrabM); 
%     end
    Rk = fun_JWR(1,0,Rs(i)); 
    Rk_r = fun_JWR(1,1,Rs(i));
    %R-D图 
    for i_wd = 1:length(wd)
        a = exp(1j*(0:Nr-1)*2*pi*fsp).';               %接收空间导向矢量
        b = exp(1j*(0:Nt-1)*2*pi*gamma*fsp).';         %发射空间导向矢量
        c = exp(1j*(0:Np-1)*2*pi*wd(i_wd)).';            %时间导向矢量
        s = kron(c,kron(b,a));                                  %遍历的导向矢量
        Pt(i,i_wd)= (abs(s'*Rk*s));                  % 匹配滤波SINR 6.27
        Ptr(i,i_wd)= (abs(s'*Rk_r*s));                  % 匹配滤波SINR 6.27
    end
end
% 
% for i_wd = 1:length(wd)
%     a = exp(1j*(0:Nr-1)*2*pi*fsp).';               %接收空间导向矢量
%     b = exp(1j*(0:Nt-1)*2*pi*gamma*fsp).';         %发射空间导向矢量
%     c = exp(1j*(0:Np-1)*2*pi*wd(i_wd)).';            %时间导向矢量
%     s = kron(c,kron(b,a));                                  %遍历的导向矢量
%     Pt(i,i_wd)= (abs(s'*Rk*s));                  % 匹配滤波SINR 6.27
%     Ptr(i,i_wd)= (abs(s'*Rk_r*s));                  % 匹配滤波SINR 6.27
% end

%% figure
% figure(1)
% [X,Y]=meshgrid(wd,RR);
% mesh(X,Y,10*log10(abs((Pt))))%10*log10
% xlabel('归一化多普勒')
% ylabel('距离/km')
% view(-0,90)
% title(['R-D图，H=',num2str(H),'m'])
figure(2)
Pt = Pt/max(max(Pt));
[X,Y]=meshgrid(wd,RR);
imagesc(wd,RR/1e3,10*log10(abs((Pt))))
xlabel('归一化多普勒')
ylabel('距离/km')
view(-0,-90)
title(['R-D图距离混叠无自转，H=',num2str(H),'m'])
figure(3)
Ptr = Ptr/max(max(Ptr));
imagesc(wd,RR/1e3,10*log10(abs((Ptr))))
xlabel('归一化多普勒')
ylabel('距离/km')
view(-0,-90)
title(['R-D图带自转距离混叠，H=',num2str(H),'m'])