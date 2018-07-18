%% 特征值
clc;clear;close all
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1250e6; %450e6   9600e6 1250e6                 %载频 Hz
C = 3e8;%299792458;                            %光速 m/s
lambda = C/fo;                       %波长 m
Nr = 16;                             %接收阵元个数
Nt = 1;                              %发射阵元个数
Np = 8;                             %相干脉冲数
D = Nt*Nr*Np;                          %自由度 
Re = 6373e3;                          %地球半径m
H = 800e3;  %700e3  %9e3                 %SBR高度 m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR速度m/s
end                
d = lambda/2;%                             %归一化阵元间距
gamma = 16;                         %发射接收空间比率
fr = 500;                         %脉冲重复频率Hz，500~2000
Tr =1/fr;                         %脉冲重复间隔s
% beta = 19.47;
beta = Vp*2/fr/d;           %%混叠系数 
alpha1 = 30;                %卫星纬度
eta = 70;                   %卫星倾角
%% 距离
Rs = 130e3;
if H>=500e3
    R =  1.149762283480509e+06;                              %地距m
    Rs = fun_R2Rs(H, R);
end

%% 距离模糊
Rsmax = fun_Rsmax(H);                       %最大探测斜距m
Rsu=C/(2*fr);                               %最大无模糊斜距m
Nk1 = floor((Rsmax-Rs)/Rsu);                %前向模糊数
Nk2 = floor((Rs-H)/Rsu);                    %后向模糊数
Rsk1 = Rs + (0:Nk1)*Rsu;                    %前向模糊斜距m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %后向模糊斜距m
Rsk = [Rsk2,Rsk1];%Rs;%                     %模糊斜距m
Rk = fun_Rs2R(H,Rsk);                       %模糊地距m
elk = fun_ELAngle(H,Rk)/180*pi;                    %各模糊距离俯仰角rad
grazek =  fun_GrazeAngle(H,Rk,Rsk)/180*pi;         %各模糊距离掠射角rad
%%
EL = asin(H/Rs);         
if H>=500e3
    EL_du = fun_ELAngle(H,R);             %波束俯仰角
    EL = EL_du/180*pi;
end

%% 偏航角幅度
CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度
%% 杂波协方差
Num = 180;
CNR = 40;
Rk = zeros(D,D);
Rk_r = zeros(D,D);
%无距离模糊
for i = 1:length(EL)%%EL
    i
    Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma,lambda);%,CrabA,CrabM 
    Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma,lambda,CrabA,CrabM);
end
E=abs(eig(Rk));
E = 10*log10(sort(E,'descend')).';
figure
hold on
plot(E,'r')
E_r=abs(eig(Rk_r));
E_r = 10*log10(sort(E_r,'descend')).';
plot(E_r,'g')
%有距离模糊
Rku = zeros(D,D);
Rku_r = zeros(D,D);
for i = 1:length(elk)%%EL
    i
    Rku_r = Rku_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(i), beta, d,gamma,lambda,CrabA,CrabM); 
    Rku = Rku + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(i), beta, d,gamma,lambda);%,CrabA,CrabM    
end
Eu=abs(eig(Rku));
Eu = 10*log10(sort(Eu,'descend')).';
hold on
plot(Eu,'b')
Eu_r=abs(eig(Rku_r));
Eu_r = 10*log10(sort(Eu_r,'descend')).';
plot(Eu_r,'k')
grid on 
legend('无距离模糊无自转','无距离模糊有自转','有距离模糊无自转','有距离模糊有自转')