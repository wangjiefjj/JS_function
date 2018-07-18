%% 全自由度stap
clc;clear;close all
isRu=0;
%% 天基平台参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1250e6; %450e6   9600e6                 %载频 Hz
c = 299792458;                            %光速 m/s
lambda = c/fo;                      %波长 m
Nr = 32;                             %接收阵元个数
Nt = 1;                              %发射阵元个数
Np = 16;                             %相干脉冲数
D = Nt*Nr*Np;                          %自由度 
Re = 6373e3;                          %地球半径m
H = 850e3;  %700e3  %9e3                 %SBR高度 m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR速度m/s
end                
d = lambda/2;%                     %阵元间距
gamma = 8;                        %发射接收空间比率
fr = 500;                         %脉冲重复频率Hz，500~2000
Tr =1/fr;                         %脉冲重复间隔s
% beta = 19.47;
beta = Vp*2/fr/d;      %%混叠系数 
alpha1 = 20/180*pi;                  %卫星纬度
eta = 90/180*pi;                    %卫星倾角

%% 偏航角幅度
CrabA = (fun_CrabAngle( alpha1,eta, H)); %偏航角
CrabM = fun_CrabMagnitude( alpha1,eta, H);%偏航幅度

%% 目标区域
R = 1000e3;                                %地距m 
Rs = fun_R2Rs(H,R);
El0 = fun_ELAngle(H,R);
% %% 模糊距离和俯仰角
Rsmax = fun_Rsmax(H);                       %最大探测斜距m
Rsu=c/(2*fr);                               %最大无模糊斜距m
%% 杂波协方差
Num = 180;                          %方位分块
CNR = 40;                          %杂噪比dB
wd = linspace(-0.5,0.5,100);       %虚拟规1化多普勒
cmj = wd./beta/(d*2*lambda);
fspc = d*cmj;
Pt = zeros(length(R),length(wd));
Ptr = zeros(length(R),length(wd));
%距离模糊
Nk1 = floor((Rsmax-Rs)/Rsu);                %前向模糊数
Nk2 = floor((Rs-H)/Rsu);                    %后向模糊数
Rsk1 = Rs + (0:Nk1)*Rsu;                    %前向模糊斜距m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %后向模糊斜距m
if isRu ==1
    Rsk = [Rsk2,Rsk1];%Rs(i);%                 %模糊斜距m
else
    Rsk = Rs;
end
RRk = fun_Rs2R(H,Rsk);                       %模糊地距m
elk = fun_ELAngle(H,RRk)/180*pi;             %各模糊距离俯仰角rad
Rk = zeros(D,D);
Rk_r = zeros(D,D);
for j = 1:length(elk)
    Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(j), beta, d,gamma, lambda);%,CrabA,CrabM
    Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(j), beta, d,gamma,lambda,CrabA,CrabM); 
end
%% Target Space-Time Steering Vector
azt = 0; elt = El0;                                      % Target azimuth and elevation angles in degrees.
fdt = 100;                                               % Target Doppler Frequency.
omegact = fdt/fr;                                        % Normalized Target Frequency.
fspt = d/lambda*cos(elt*pi/180)*cos(azt*pi/180);     % Target Spatial Frequency.
at = exp(1i*2*pi*fspt*(0:Nr-1));                          % Target Spatial Steering Vector.
bt = exp(1i*2*pi*omegact*(0:Np-1));                       % Target Doppler Steering Vector
vt = kron(bt,at).';                                      % Target Space-Time Steering Vector.

%% Optimum, Fully Adaptive STAP Solution
w = Rk\vt;                                               % Eq. (104)
w_r = Rk_r\vt;
%% Adapted Patterns
az = -90:.5:90;     Laz = length(az);
fd = -150:.5:150;  Lfd = length(fd);
fsp = d/lambda*cos(El0)*cos(az*pi/180);
omega = fd/fr;
Pw1 = zeros(Lfd,Laz);
for m=1:Laz
    for n=1:Lfd
        a = exp(1i*2*pi*fsp(m)*(0:Nr-1));                % Dummy Spatial Steering Vector.(Dummy虚拟)
        b = exp(1i*2*pi*omega(n)*(0:Np-1));              % Dummy Doppler Steering Vector
        v = kron(b,a).';
        Pw1(n,m) = abs(w'*v)^2;
    end
end

%% Normalization:
max_value = max(max(Pw1));
Pw = Pw1/max_value;
%% Plot the Adapted Pattern
figure('NumberTitle', 'off','Name', ...
       'Figure 23a. Example Scenario: Adapted Pattern for Optimum Fully Adaptive STAP', ...
       'Position',[1 1 700 600]);
[Az Doppler] = meshgrid(sin(az*pi/180),fd);
colormap jet;
mesh(Az, Doppler, 10*log10(abs(Pw)));
shading interp;
xlim([-1 1])
ylim([-150 150]);
xlabel('sin(Azimuth)');
ylabel('Doppler Frequency (Hz)');
h = colorbar;
% h = colorbar('YTickLabel',{-80:10:0});
set(get(h,'YLabel'),'String','Relative Power (dB)');