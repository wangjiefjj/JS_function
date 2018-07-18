%% 根据 天基雷达空时自适应杂波抑制技术,张增辉 参数仿的
clc
clear
close all
% isRu: 有无距离模糊，有1，无0
% isRotation：有无自转，有1，无0
isRu = 0;
isRotation = 0;
%% 常数
Re = 6373e3;                    %圆地球半径m
C  = 3e8;                 %光速m/s 299792458
%% 雷达系统参数
fo = 1250e6;                    %载频Hz
lambda = C/fo;                  %波长
Nr_az = 16;%16                     %接收方位子阵个数
Nr_el = 4;%25                     %接收俯仰子阵个数
Nt_az = 16;%16                     %发射方位子阵个数
Nt_el = 4;%25                     %发射俯仰子阵个数
Ti = 14.2e-3;                   %驻留时间s
fr = 500;                      %脉冲重复周期Hz
Np = 16;%floor(fr*Ti);              %接收脉冲数 
Pt = 30e3;                      %峰值功率W
Tp = 20e-6;                     %脉宽s
Pa = Pt*Tp/(1/fr);              %平均发射功率
Ls = 6;                         %系统损耗dB
B = 1e6;                        %带宽Hz
D = Tp*B;                       %脉压比
Gt = 40;                        %发射增益dB
Gr = 23;                        %接收增益dB
d = lambda/2;                   %阵元间距m
d_saz = 0.12*25;                %方位向子阵间距m
L = 48;                         %天线长度m
W = 3;                          %天线宽度m
%% 雷达
% alpha_r = 30;                   %雷达纬度N,deg
% beta_r = 110;                   %雷达经度E,deg
H = 850e3;                      %高度m
Vp = fun_Vp(H);                 %平台速度
eta = 70;                       %轨道倾角,deg
alpha1 = 110;                   %雷达经度E,deg
beta1 =  30;                    %雷达纬度N,deg
beta = Vp*2/fr/d;               %%混叠系数;
%% 目标
graze = 30;                     %天线法线掠射角deg
alpha2 = 120;                   %目标经度E,deg
beta2  = 26;                    %目标纬度N,deg
RCS = 5;                        %目标RCS，dBsm，（dB*m^2）
% 根据graze得到的地距，即法线方向的波束指向地距
R_graze = ( pi/2 - graze/180*pi - asin( Re/(Re+H) * cos(graze/180*pi) ) ) * Re; %(2.26)变形
Rs_graze = fun_R2Rs(H,R_graze);
% 方位角
Nc = 180;                           %杂波块个数 
Az = linspace(0,179,Nc);            %方位角度
El = linspace(0,179,Nc);            %俯仰向角度
% Az = Az/180*pi;
% El = El/180*pi;
EL =  fun_ELAngle(H,R_graze);       %倾斜角（天线法线俯仰角）
LAz = length(Az);
AFt = ones(1,LAz);                  % 发射天线阵因子
AFr = ones(1,LAz);                  % 接收天线阵因子
%% 波束指向计算
All0 = fun_ComputeAzEl(alpha1,beta1,H,eta,alpha2,beta2,0);

[X1,Y1,Z1] = fun_JWH2XYZ(alpha1,beta1, H+Re); %卫星经纬高转XYZ
[X2,Y2,Z2] = fun_JWH2XYZ(alpha2,beta2, Re);   %目标点经纬高转XYZ
Rs = sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);     %目标斜距
R = fun_Rs2R(H,Rs);                           %目标地距
% %计算波束指向的俯仰角
% graze0 = fun_GrazeAngle(H,R,Rs);
% El0 = fun_ELAngle(H,R);
% % 计算波束指向的方位角
% % 卫星雷达坐标系坐标
% mu = asin(sin(beta1/180*pi)/sin(eta/180*pi));
% [Xr1,Yr1,Zr1] = fun_XYZ2Radar(alpha1,beta1,eta,Re+H, X1,Y1,Z1);
% Az0 = abs(acos(Xr1/Rs/sin(El0/180*pi)))/pi*180;
%% 方向图 (看看方向图啥样，后面计算还要重新算)
Iraz = taylorwin(Nr_az).';                   %%方位接收泰勒加权
Irel = taylorwin(Nr_el).';                   %%俯仰接收泰勒加权
Itaz = taylorwin(Nt_az).'; %,4,-23         %%方位发射泰勒加权
Itel = taylorwin(Nt_el).'; %,4,-23         %%俯仰发射泰勒加权
% It1 = ones(1,Nr_az);                        %%发射均匀加权
% It2 = ones(1,Nr_el);                        %%发射均匀加权
Az0 = All0.Az;   %91.42                       %% 主波束方向
El0 = All0.El;   %50.18                       %%俯仰角
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%掠射角(2.27)deg
% graze0 = acos(sin(El0/180*pi)*(Re+H)/Re)/pi*180; 
%(2.26)变形
% R = ( pi/2 - graze0/180*pi - asin( Re/(Re+H) * cos(graze0/180*pi) ) ) * Re; 
% Rs = fun_R2Rs(H,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fr = zeros(length(Az),length(El));          %接收方向图
Ft = zeros(length(Az),length(El));          %发射方向图
% % 每一列是一个方位角，每一行是一个俯仰角
for i=1:length(Az)
%     i 
    for j = 1:length(El)% 
        Fr(j,i) = fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Iraz,Irel);
%         Fr(j,i) = fun_Ax(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,Ir1,Ir2);
        Ft(j,i) = fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Itaz,Itel);
%         Ft(j,i) = fun_Ax(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,It1,It2);
    end
end
Fr = 10*log10(abs(Fr));
Fr(Fr<-40)=-40;
Ft = 10*log10(abs(Ft));
Ft(Ft<-40)=-40;
% %%接收，画图
figure()
mesh(Az,El,Fr)%10*log10
xlabel('方位角/deg')
ylabel('俯仰角/deg')
zlabel('幅度/dB')
title('接收方向图')
% hold on
% plot(Az0,El0,'r.','MarkerSize',20)
% [maxAz_r,maxEl_r] = find(Fr == max(max(Fr)));
% figure()
% plot(Az,Fr(:,maxEl_r))
% title('接收俯仰向方向图')
% xlabel('幅度/deg')
% ylabel('方位角/dB')
% figure()
% plot(El,Fr(maxAz_r,:))
% title('接收方位方向图')
% xlabel('幅度/dB')
% ylabel('俯仰角/deg')
% % % 发射，画图
% [label_az, label_el] = meshgrid(El,Az);
% figure()
% mesh(label_az,label_el,Ft)%10*log10
% xlabel('方位角/deg')
% ylabel('俯仰角/deg')
% zlabel('幅度/dB')
% title('发射方向图')
% figure()
% imagesc(El,Az,Ft)%10*log10
% xlabel('方位角/deg')
% ylabel('俯仰角/deg')
% title('发射方向图')
% [maxEl_t,maxAz_t] = find(Ft == max(max(Ft)));
% figure()
% plot(Az,Ft(maxEl_t,:))
% title('发射方位向方向图')
% xlabel('方位角/deg')
% ylabel('幅度/dB')
% figure()
% plot(El,Ft(:,maxAz_t))
% hold on 
% plot(El,Ft(:,maxAz_r),'r')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 根据目标经纬度画出波束覆盖威力图
% LL = 200;
% alpha2_t = linspace(alpha2-20,alpha2+20,LL);
% beta2_t = linspace(beta2-20,beta2+20,LL);
% Fr = zeros(LL,LL);
% Ft = zeros(LL,LL);
% for i_a = 1:LL
%     i_a
%     for i_b = 1:LL
%         All(i_a,i_b) = fun_ComputeAzEl(alpha1,beta1,H,eta,alpha2_t(i_a),beta2_t(i_b),0);
% %         Fr(i_a,i_b) = fun_F(Nr_az,Nr_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
% %             All0.Az,All0.El,EL,Ir1,Ir2);
%         Fr(i_a,i_b) = fun_Ax(Nr_az,Nr_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
%             All0.Az,All0.El,Ir1,Ir2);
% %         Ft(i_a,i_b) = fun_F(Nt_az,Nt_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
% %             All0.Az,All0.El,EL,It1,It2);
%         Ft(i_a,i_b) = fun_Ax(Nt_az,Nt_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
%             All0.Az,All0.El,It1,It2);
%     end
% end
% 
% 
% Fr = 10*log10(abs(Fr));
% % Fr(Fr<-40)=-40;
% Ft = 10*log10(abs(Ft));
% Ft(Ft<-40)=-40;
% 
% figure
% imagesc(beta2_t,alpha2_t,abs(Fr))
% xlabel('纬度')
% ylabel('经度')
% hold on 
% plot(beta1,alpha1,'g.','MarkerSize',20) %%雷达位置
% hold on 
% plot(beta2,alpha2,'r.','MarkerSize',20) %%主波束指向 
% % axis([beta2-5,beta2+5,alpha2-10,alpha2+5,-20,25])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 终端SNR计算
k = 1.3806488e-23;                  % 玻尔兹曼常数 J/K.
T0 = 300;                           % 标准开尔文室温 Kelvin.
Fn = 3;                             % 接收机噪声系数 dB;
Te = T0*(10^(Fn/10));               % Effective Receiver Temperature in Kelvin.
Ts = 10^(Ls/10)*Te;                 % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                          % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                          % Receiver Noise Power in Watts
sigma2 = 1;                         % Normalized Noise Power in Watts.
S = 350e3*350e3;                    %探测面积 (自定义)
dS = Rs*lambda/L*Rs*lambda/W;       %每个波位的面积 根据P17
A = L*W;                           %天线孔径m
SNR = ( Pa*A^2*10^(RCS/10)*Ti ) / ( 4*pi*lambda^2*k*T0*10^(Fn/10)*Rs^4*10^(Ls/10));%(2.21)
SNR = 10*log10(SNR);

%% 功率孔径积
PA = 4*pi*k*T0*Fn*Rs^2*10^(SNR/10)*S*sin(graze/180*pi)/10^(RCS/10);
PA = 10*log10(PA);

%% 距离模糊
Rsmax = fun_Rsmax(H);                       %最大探测斜距m
r = C*Tp/2/D;                               %距离分辨率m
Rsu=C/(2*fr);                               %最大无模糊斜距m
Nk1 = floor((Rsmax-Rs)/Rsu);                %前向模糊数
Nk2 = floor((Rs-H)/Rsu);                    %后向模糊数
Rsk1 = Rs + (0:Nk1)*Rsu;                    %前向模糊斜距m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %后向模糊斜距m
if isRu==1
    Rsk = [Rsk2,Rsk1];%Rs;%                 %模糊斜距m
elseif isRu == 0
    Rsk = Rs;%
end
Rk = fun_Rs2R(H,Rsk);                       %模糊地距m
elk = fun_ELAngle(H,Rk);                    %各模糊距离俯仰角deg
grazek =  fun_GrazeAngle(H,Rk,Rsk);         %各模糊距离掠射角deg
% grazek = acos(sin(elk/180*pi)*(Re+H)/Re)/pi*180; 
sk = r./cos(grazek/180*pi)*2*pi.*Rk/Nc;     % 杂波单元面积m^2  (2.36)
% sak = sk*Nc;                              %杂波环面积m^2
%% 后项散射系数
%<<地杂波背景中机载预警雷达作用距离分析>>
for i = 1:length(grazek)
    sigma0(1,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,1));     %沙漠
    sigma0(2,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,2));     %农田
    sigma0(3,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,3));     %丘陵
    sigma0(4,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,4));     %高山
    sigma0(5,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,5,2));   %海洋
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
%% 各杂波块雷达等效散射面积 
sigmac = 10.^(sigma0(1,:)./10).*sk;                     %雷达等效散射面积(2.35)
sigmac = 10*log10(sigmac);                              %dB
%% 每个杂波单元的CNR
%  计算接收和发射增益
for i = 1:Nc
    for j = 1:length(elk)
        Ft_t=abs(fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,Itaz,Itel)).^2;
%         Ft_t=abs(fun_Ax(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,Itaz,Itel)).^2;
        Gtin(j,i) = Gt*Ft_t; 
        Fr_t=abs(fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,Iraz,Irel)).^2;
%         Fr_t=abs(fun_Ax(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,Iraz,Irel)).^2;
        Grin(j,i) = Gr*Fr_t;
    end
end

% 计算每个方位角(杂波块)的CNR
for i = 1:Nc
    for j = 1:length(elk)
         CNR(j,i) = (Pt*Gtin(j,i)*Grin(j,i)*lambda^2*10.^(sigmac(j)/10)/((4*pi)^3*Pn*10^(Ls/10)*(Rsk(j)^4)'));
    end   
end
CNR = sigma2 * CNR; %
% CNR = 10 * log10(CNR);
index = find(CNR==max(CNR));
% CNR(index)=1000000;
% CNR(1:index-1) = 1e-4;
% CNR(index+1:end) = 1e-4;
% figure
% mesh((CNR))
% figure()
% [X,Y] = meshgrid(Az,elk);
% mesh(X,Y,CNR)
% xlabel('方位角/deg')
% ylabel('俯仰角/deg')
% zlabel('CNR/dB')

%% 杂波协方差计算，计算每个杂波块的空间和时间频率（elk，azk）
azk = Az;
cmj1 = sin(elk/180*pi).' * cos(azk/180*pi);                    %入射锥角
fspc1 = d/lambda*cmj1;                                       %空间频率
cmj2 = sin(elk/180*pi).' * sin(azk/180*pi);                    %入射锥角
fspc2 = d/lambda*cmj2;                                       %空间频率
% 航偏角，航偏幅度
if isRotation == 1
    CrabA = fun_CrabAngle( beta1,eta, H);%0;%                      %偏航角/rad
    CrabM = fun_CrabMagnitude( beta1,eta, H); %1;%                  %偏航幅度
elseif isRotation == 0
    CrabA = 0;%                  %偏航角/rad
    CrabM = 1;%                  %偏航幅度
end
omegac = beta * d/lambda*2 * CrabM * (sin(elk/180*pi).'*cos(azk/180*pi+CrabA));   %杂波归1化多普勒
% mod(omegac(index),0.5)
Vc = zeros(Nr_az*Np*Nr_el,Nc);
% Vc = zeros(Nr_az*Np,Nc);
Rc = zeros(Nr_az*Np*Nr_el,Nr_az*Np*Nr_el);
% load Ksic.mat
for j = 1:length(elk)
    j
    for k=1:Nc
        % 空间导向矢量.
        a(:,k) = exp(1i*2*pi*fspc1(j,k)*(0:Nr_az-1));    % 空间导向矢量,方位向/sqrt(Nr_az)
        b(:,k) = exp(1i*2*pi*omegac(j,k)*(0:Np-1));     % Temporal Steering Vector/sqrt(Np)
        c(:,k) = exp(1i*2*pi*fspc2(j,k)*(0:Nr_el-1));  % 空间导向矢量,俯仰向/sqrt(Nr_el)
        Vc(:,k) = kron(b(:,k),kron(c(:,k),a(:,k))); 
    end
    
    A = diag(fliplr(CNR(j,:)));
%     A = repmat(CNR(j,:),Nr_az*Np*Nr_el,1);
    Rc = Rc + Vc*A*Vc';
end
 
Rn = sigma2*eye(Nr_az*Np*Nr_el);
% Rn = sigma2*eye(Nr_az*Np);
Rcn = Rc + Rn;
E=abs(eig(Rcn));
E = 10*log10(sort(E,'descend')).';
figure
hold on
plot(E(1:100),'r')
figure(3)
mesh(abs(Rcn))
%% MVD-plot
L2 = 500;
%%归一化时间频率 
omega = linspace(-0.5,0.5,L2);
fsp1 = 0;%d/lambda/2*cmj; 
fsp2 = d/lambda*sin(El0/180*pi);%d/lambda/2*cmj; 
iRcn = inv(Rcn);
for j = 1:length(omega)  
    a = exp(1i*2*pi*fsp1*(0:Nr_az-1)).';    % 空间导向矢量./sqrt(Nr_az)
    c = exp(1i*2*pi*fsp2*(0:Nr_el-1)).';%/sqrt(Nr_el)
    b = exp(1i*2*pi*omega(j)*(0:Np-1)).';    % 时间导向矢量/sqrt(Np)
%      v = kron(b,a);                           % 空时导向矢量.
     v = kron(b,kron(c,a));
     MVD(j) = abs(v'*iRcn*v)^2;
end
figure()
plot(omega,10*log10(abs(MVD)/max(max(abs(MVD)))))%10*log10(abs(MVD)/max(max(abs(MVD))))
xlabel('归一化时间频率')
ylabel('SINR/dB')
% axis([-0.5,0.5,-60,-40])
%% MVD——mesh
% L2 = 200;
% fsp = linspace(-0.5,0.5,L2);
% omega = linspace(-0.5,0.5,L2);              %% 归一化时间频率 
% iRcn = inv(Rcn);
% for i = 1:L2
%     i
%     for j = 1:L2
%         a = exp(1i*2*pi*fsp(i)*(0:Nr_az-1)).';    % 空间导向矢量.
%         b = exp(1i*2*pi*omega(j)*(0:Np-1)).';    % 时间导向矢量
%         v = kron(b,a);                           % 空时导向矢量.
%         MVD(i,j) = abs(v'*iRcn*v)^2;
%     end
% end
% 
% figure()
% [X,Y] = meshgrid(omega,fsp);
% mesh(X,Y,-10*log10(abs(MVD)/max(max(abs(MVD)))))
% xlabel('归一化空间频率')
% ylabel('归一化时间频率')
