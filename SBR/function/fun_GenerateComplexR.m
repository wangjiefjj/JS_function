function [ Rcn,El0,d,lambda,fr] = fun_GenerateComplexR( isRu,isRotation )
%   产生需要的杂波协方差，通过详细数据来的
% isRu: 有无距离模糊，有1，无0
% isRotation：有无自转，有1，无0
%% 常数
Re = 6373e3;                    %圆地球半径m
C  = 3e8;                       %光速m/s 299792458
%% 雷达系统参数
fo = 1.25e9;                    %载频Hz
lambda = C/fo;                  %波长
Nr_az = 8;%16                     %接收方位子阵个数
Nr_el = 8;%25                     %接收俯仰子阵个数
Nt_az = 8;%16                     %发射方位子阵个数
Nt_el = 8;%25                     %发射俯仰子阵个数
Ti = 14.2e-3;                   %驻留时间s
fr = 10000;                      %脉冲重复周期Hz
Np = 8;%floor(fr*Ti);              %接收脉冲数 
Pt = 30e3;                      %峰值功率W
Tp = 20e-6;                     %脉宽s
Pa = Pt*Tp/(1/fr);              %平均发射功率
Ls = 6;                         %系统损耗dB
B = 1e6;                        %带宽Hz
D = Tp*B;                       %脉压比
Gt = 23;                        %发射增益dB
Gr = 23;                        %接收增益dB
d = lambda/2;                   %阵元间距m(俯仰向)
% d_saz = d*25;                %方位向子阵间距m
L = 48;                         %天线长度m
W = 3;                          %天线宽度m
%% 雷达
% alpha_r = 30;                   %雷达纬度N,deg
% beta_r = 110;                   %雷达经度E,deg
H = 850e3;                        %高度m
Vp = fun_Vp(H);                   %平台速度
eta = 70;                         %轨道倾角,deg
alpha1 = 110;                     %雷达经度E,deg
beta1 =  30;                      %雷达纬度N,deg
beta = Vp*2/fr/d;                 %%混叠系数;
%% 目标
graze = 30;                     %掠射角deg
alpha2 = 120;                   %目标经度E,deg
beta2  = 25;                    %目标纬度N,deg
RCS = 5;                        %目标RCS，dBsm，（dB*m^2）
% 根据graze得到的地距，即法线方向的波束指向直到目标时的地距
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

%% 设定方向图 (看看方向图啥样，后面计算还要重新算)
Az0 = 91.671396591855710;   %91.42                       %% 主波束方向
El0 = 49.836208625040160;   %50.18                       %%俯仰角
%%掠射角(2.27)deg
graze0 = acos(sin(El0/180*pi)*(Re+H)/Re)/pi*180; 
%(2.26)变形
R = ( pi/2 - graze0/180*pi - asin( Re/(Re+H) * cos(graze0/180*pi) ) ) * Re; 
Rs = fun_R2Rs(H,R);
% 方向图计算
% Ir1 = taylorwin(Nr_az).';                   %%方位接收泰勒加权
% Ir2 = taylorwin(Nr_el).';                   %%俯仰接收泰勒加权
Ir1 = ones(1,Nr_az);                      %%接收均匀加权
Ir2 = ones(1,Nr_el);                      %%接收均匀加权
% It1 = taylorwin(Nr_az,4,-23).'; %,4,-23           %%方位发射泰勒加权
% It2 = taylorwin(Nr_el,4,-23).'; %,4,-23           %%俯仰发射泰勒加权
It1 = ones(1,Nr_az);                      %%发射均匀加权
It2 = ones(1,Nr_el);                      %%发射均匀加权
Fr = zeros(length(Az),length(El));          %接收方向图
Ft = zeros(length(Az),length(El));          %发射方向图
% 每一列是一个方位角，每一行是一个俯仰角
% for i=1:length(Az)
% %     i
%     for j = 1:length(El)% 
%         Fr(j,i) = fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Ir1,Ir2);
%         Ft(j,i) = fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,It1,It2);
%     end
% end
% Fr = 10*log10(abs(Fr*25));
% Fr(Fr<-40)=-40;
% Ft = 10*log10(abs(Ft*25));
% Ft(Ft<-40)=-40;
% % %%接收，画图
% % [label_az, label_el] = meshgrid(El,Az);
% % figure()
% % mesh(label_az,label_el,Fr)%10*log10
% % xlabel('方位角/deg')
% % ylabel('俯仰角/deg')
% % zlabel('幅度/dB')
% % title('接收方向图')
% [maxAz_r,maxEl_r] = find(Fr == max(max(Fr)));
% figure()
% plot(Az,Fr(:,maxEl_r))
% title('接收方位向方向图')
% xlabel('幅度/deg')
% ylabel('方位角/dB')
% figure()
% plot(El,Fr(maxAz_r,:))
% title('接收俯仰向方向图')
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
    disp('isRu')
    Rsk = [Rsk2,Rsk1];%Rs;%                 %模糊斜距m
elseif isRu == 0
    disp('noRu')
    Rsk = Rs;%
end
Rk = fun_Rs2R(H,Rsk);                       %模糊地距m
elk = fun_ELAngle(H,Rk);                    %各模糊距离俯仰角deg
grazek =  fun_GrazeAngle(H,Rk,Rsk);         %各模糊距离掠射角deg
% grazek = acos(sin(elk/180*pi)*(Re+H)/Re)/pi*180; 
sk = r./cos(grazek/180*pi)*2*pi.*Rk/Nc*10;     % 杂波单元面积m^2  (2.36)
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
        Ft_t=abs(fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,It1,It2)).^2;
        Gtin(j,i) = Gt*Ft_t; 
        Fr_t=abs(fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,Ir1,Ir2)).^2;
        Grin(j,i) = Gr*Fr_t;
    end
end
% 计算每个方位角(杂波块)的CNR
for i = 1:Nc
    for j = 1:length(elk)
%          CNR(j,i) = sum(Pt.*Gtin(j,i).*Grin(j,i)*lambda^2.*10.^(sigmac.'/10)./((4*pi)^3*Pn*10^(Ls/10).*(Rsk(j).^4).'));
          CNR(j,i) = (Pt*Gtin(j,i)*Grin(j,i)*lambda^2*10.^(sigmac(j)/10)/((4*pi)^3*Pn*10^(Ls/10)*(Rsk(j)^4)'));
    end   
end
CNR = sigma2 * CNR; %
CNR = 10 * log10(CNR);

% figure()
% [X,Y] = meshgrid(Az,elk);
% mesh(X,Y,CNR)
% xlabel('方位角/deg')
% ylabel('俯仰角/deg')
% zlabel('CNR/dB')

%% 杂波协方差计算，计算每个杂波块的空间和时间频率（elk，azk）
azk = Az;
% 航偏角，航偏幅度
if isRotation == 1
    disp('isRotation')
    CrabA = fun_CrabAngle( beta1,eta, H);%0;%                      %偏航角/rad
    CrabM = fun_CrabMagnitude( beta1,eta, H); %1;%                  %偏航幅度
elseif isRotation == 0
    disp('noRotation')
    CrabA = 0;%                  %偏航角/rad
    CrabM = 1;%                  %偏航幅度
end


Vc = zeros(Nr_az*Np*Nr_el,Nc);
% Vc = zeros(Nr_az*Np,Nc);
Rc = zeros(Nr_az*Np*Nr_el,Nr_az*Np*Nr_el);
for j = 1:length(elk)
    j
    for k=1:Nc
        cmj = cos(elk(j)/180*pi).' * cos(azk(k)/180*pi);                    %入射锥角
        fspc = d/lambda*cmj;                                       %空间频率
        omegac = beta * d/lambda*2 * CrabM * (cos(elk(j)/180*pi).*cos(azk(k)/180*pi+CrabA));   %杂波归1化多普勒
        % 空间导向矢量.
        a = exp(1i*2*pi*fspc*(0:Nr_az-1)).';    % 空间导向矢量,方位向
        c = exp(1i*2*pi*fspc*(0:Nr_el-1)).';    % 空间导向矢量,俯仰向
        b = exp(1i*2*pi*omegac*(0:Np-1)).';     % Temporal Steering Vector 
        Vc(:,k) = kron(b,kron(c,a));
    end
    A = diag((CNR(j,:)));
    Rc = Rc+Vc*A*Vc';
end
% Rc = Nr_az*Np*Nr_el*Rc/sum(eig(Rc));
Rn = sigma2*eye(Nr_az*Np*Nr_el);
% Rn = sigma2*eye(Nr_az*Np);

Rcn = Rc + Rn;


end

