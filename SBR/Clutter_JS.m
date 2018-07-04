%% 根据 天基雷达空时自适应杂波抑制技术,张增辉 参数仿的
clc
clear
close all
%% 常数
Re = 6373e3;                    %圆地球半径m
c  = 3e8;                 %光速m/s 299792458
%% 雷达系统参数
fo = 1.25e9;                    %载频Hz
lambda = c/fo;                  %波长
Nr_az = 16;                     %接收方位子阵个数
Nr_el = 25;                     %接收方位子阵个数
Nt_az = 16;                     %发射方位子阵个数
Nt_el = 25;                     %发射方位子阵个数
Ti = 14.2e-3;                   %驻留时间s
fr = 5000;                      %脉冲重复周期Hz
Np = fr*Ti;                     %接收脉冲数 
Pt = 30e3;                      %峰值功率W
Tp = 20e-6;                     %脉宽s
Pa = Pt*Tp/(1/fr);              %平均发射功率
Ls = 6;                         %系统损耗dB
B = 1e6;                        %带宽Hz
D = Tp*B;                       %脉压比
Gt = 23;                        %发射增益dB
Gr = 23;                        %接收增益dB
d = lambda/2;                   %阵元间距m(俯仰向)
d_saz = 0.12*25;                %方位向子阵间距m
%% 雷达
alpha_r = 30;                   %雷达纬度N,deg
beta_r = 110;                   %雷达经度E,deg
H = 850e3;                      %高度m
Vp = fun_Vp(H);                 %平台速度
phi = 70;                       %轨道倾角 
alpha1 = 110;                   %雷达经度E,deg
beta1 =  30;                    %雷达纬度N,deg
%% 目标
graze = 30;                     %掠射角deg
alpha2 = 120;                  %目标经度E,deg
beta2  = 25;                   %目标纬度N,deg
RCS = 5;                        %目标RCS，dBsm，（dB*m^2）
R = ( pi/2 - graze/180*pi - asin( Re/(Re+H) * cos(graze/180*pi) ) ) * Re; %(2.26)变形
Rs = fun_R2Rs(H,R);
% 方位角
Nc = 359;                           %杂波块个数 
Az = linspace(0,179,Nc);            %方位角度
El = linspace(0,179,Nc);            %俯仰向角度
% Az = Az/180*pi;
% El = El/180*pi;
EL =  fun_ELAngle(H,R);              %倾斜角（天线法线俯仰角）
LAz = length(Az);
AFt = ones(1,LAz);                  % 发射天线阵因子
AFr = ones(1,LAz);                  % 接收天线阵因子

%% 终端SNR计算
k = 1.3806488e-23;                  % 玻尔兹曼常数 J/K.
T0 = 300;                           % 标准开尔文室温 Kelvin.
Fn = 3;                             % 接收机噪声系数 dB;
S = 350e3*350e3;                    %探测面积 (自定义)
dS = Rs*lambda/48*Rs*lambda/3;      %每个波位的面积

A = 3*48;                           %天线孔径m
SNR = ( Pa*A^2*10^(RCS/10)*Ti ) / ( 4*pi*lambda^2*k*T0*10^(Fn/10)*Rs^4*10^(Ls/10) );
SNR = 10*log10(SNR);

%% 功率孔径积
PA = 4*pi*k*T0*Fn*Rs^2*10^(SNR/10)*S*sin(graze/180*pi)/10^(RCS/10);
PA = 10*log10(PA);

%% 方向图 
Az0 = 91.42;                 %% 主波束方向
El0 = 50.18;
Ir1 = taylorwin(Nr_az).';    %%接收泰勒加权
Ir2 = taylorwin(Nr_el).';    %%接收泰勒加权
It1 = taylorwin(Nr_az,4,-23).';    %%发射泰勒加权
It2 = taylorwin(Nr_el,4,-23).';    %%发射泰勒加权
% It1 = ones(1,Nr_az);    %%发射均匀加权
% It2 = ones(1,Nr_el);    %%发射均匀加权
Fr = zeros(length(Az),length(El));  %接收方向图
Ft = zeros(length(Az),length(El));  %发射方向图
for i=1:length(Az)
    i
    for j = 1:length(El)% 
        Fr(i,j) = fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Ir1,Ir2);
        Ft(i,j) = fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,It1,It2);
    end
end
Fr = 10*log10(abs(Fr*25));
Fr(Fr<-40)=-40;
Ft = 10*log10(abs(Ft*25));
Ft(Ft<-40)=-40;
% % 接收，画图
% [label_az, label_el] = meshgrid(El,Az);
% figure()
% mesh(label_az,label_el,Fr)%10*log10
% xlabel('俯仰角/deg')
% ylabel('方位角/deg')
% zlabel('幅度/dB')
% title('接收方向图')
% [maxAz_r,maxEl_r] = find(Fr == max(max(Fr)));
% figure()
% plot(Az,Fr(:,maxEl_r))
% title('接收方位向方向图')
% xlabel('方位角/deg')
% ylabel('幅度/dB')
% figure()
% plot(El,Fr(maxAz_r,:))
% title('接收俯仰向方向图')
% xlabel('俯仰角/deg')
% ylabel('幅度/dB')
% % % 发射，画图
% [label_az, label_el] = meshgrid(El,Az);
% figure()
% mesh(label_az,label_el,Ft)%10*log10
% xlabel('俯仰角/deg')
% ylabel('方位角/deg')
% zlabel('幅度/dB')
% title('发射方向图')
% figure()
% imagesc(El,Az,Ft)%10*log10
% xlabel('俯仰角/deg')
% ylabel('方位角/deg')
% title('发射方向图')
% [maxAz_t,maxEl_t] = find(Ft == max(max(Ft)));
% figure()
% plot(Az,Ft(:,maxEl_t))
% title('发射方位向方向图')
% xlabel('方位角/deg')
% ylabel('幅度/dB')
% figure()
% plot(El,Ft(maxAz_t,:))
% hold on 
% plot(El,Ft(maxAz_r,:),'r')
% title('发射俯仰向方向图')
% xlabel('俯仰角/deg')
% ylabel('幅度/dB')
% figure()
% index = find(Ft<10);
% label_az_t = label_az;
% label_az_t(index)=nan;
% label_el_t = label_el;
% label_el_t(index)=nan;
% Ft_t = Ft;
% Ft_t(index)=nan;
% contour(label_az_t,label_el_t,Ft_t)
%% 距离模糊
Rsmax = fun_Rsmax(H);                   %最大探测斜距m
r = c*Tp/2/D;                           %距离分辨率m
Rsu=c/(2*fr);                           %最大无模糊斜距m
Nk1 = floor((Rsmax-Rs)/Rsu);            %前向模糊数
Nk2 = floor((Rs-H)/Rsu);                %后向模糊数
Rsk1 = Rs + (0:Nk1)*Rsu;                %前向模糊斜距m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;             %后向模糊斜距m
Rsk = [Rsk2,Rsk1];                      %模糊斜距m
Rk = fun_Rs2R(H,Rsk);                   %模糊地距m
grazek =  fun_GrazeAngle(H,Rk,Rsk);     %各模糊距离掠射角deg
sk = r./cos(grazek/180*pi)*2*pi.*Rk/Nc; % 杂波块面积m^2
% sak = sk*Nc;                            %杂波环面积m^2
%% 后项散射系数
%<<地杂波背景中机载预警雷达作用距离分析>>
for i = 1:length(grazek)
    sigma0(1,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,1));
    sigma0(2,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,2));
    sigma0(3,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,3));
    sigma0(4,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,4));  
    sigma0(5,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,5,2));
end
% figure()
% hold on
% plot(grazek,sigma0(1,:),'r')
% plot(grazek,sigma0(2,:),'g')
% plot(grazek,sigma0(3,:),'b')
% plot(grazek,sigma0(4,:),'k')
% plot(grazek,sigma0(5,:),'c')
% grid on
% box on
% legend('沙漠','农田','丘陵','高山','海洋')
% axis([0,76,-80,20]);
% xlabel('掠射角/deg')
% ylabel('RCS(\sigma_0)/dB')
% title('杂波后向散射系数随掠射角的变化')

