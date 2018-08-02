function [ P ] = fun_GenerateComplexTrainData( isRu,isRotation,Num_TrainData)
%% 产生不同距离单元的杂波数据以及杂波协方差

%% 本程序用于产生机载雷达回波仿真数据
%% 生成仿真参数
%% 雷达参数
if nargin == 2
    Num_TrainData = 1;
end
N=8;        %每行阵元数
M=8;        %每列阵元数
K=8;        %脉冲数
lambda=0.5; %波长
c=3e8;      %光速
d=0.25;     %阵元间距
fr=10000;   %脉冲重复频率
% Rmax=para_data(10);%最大作用距离 【无用参数，后面根据场景计算】

Inn=30;             %行向所加切比雪夫权重dB
Imm=30;             %列向所加切比雪夫权重dB
Ikk=30;             %时域所加切比雪夫权重dB
% CNR_dB=para_data(16);%杂噪比dB
% CNR=10^(CNR_dB/10);%杂噪比
thetap=0;           %阵面和航迹之间夹角
senser_error=0;     %阵元误差
channel_error=0;    %通道误差
%channel_error=0.05;%通道误差
gamma=1;
xinhao=1;
Q=1;                %% 信号数量,MIMO的话Q=8
%% 平台参数
H=800e3;                    %高度
V=fun_Vp(H);                %速度
rsmax = fun_Rsmax(H);
rs=H;%1300e3;          %
r=fun_Rs2R(H,rs);
el0=fun_ELAngle(H,r)/180*pi;        %波束指向俯仰角
az0=pi/2;                           %波束指向方位角
alpha1 = 30;                        %平台纬度
eta = 70;                           %平台倾角
if isRotation==1
    crabA = fun_CrabAngle(alpha1,eta, H);       %偏航角rad
    crabM = fun_CrabMagnitude(alpha1,eta, H);   %%偏航幅度
else
    crabA = 0;      %偏航角rad
    crabM = 1;      %%偏航幅度 
end

%% 杂波
Br=0;               %杂波起伏
flag_clutter=5;     %杂波类型，分别为海杂波、沙漠、农田、丘陵、高山
flag_sea=5;         %海情等级1-5
flag_kind=3;        %杂波幅度分布类型，分别为Log-normal、Weibull和K分布
flag_noise=0;       %是否含有噪声
Rs=eye(M);
% Rs=eye(M)+tril(0.8*rand(M,M),-1)+triu(0.8*rand(M,M),1);
% Rs=eye(M)+tril(0.05*rand(M,M),-1)+triu(0.05*rand(M,M),1);
%====设置参数
N_theta = 180; % ####方位角所分份数
N_psi=180;      %锥角所分份数&& 没用到的设置
dpsi=pi/180;    %锥角步长&& 没用到的设置
N_phi=90;       %俯仰角所分份数&& 没用到的设置
dphi=pi/90;     %俯仰角步长&& 没用到的设置
dtheta=pi/180;  %方位角步长&& 没用到的设置
% % caa=(channel_error*randn(N,1));ca=caa*caa.';%通道间幅度误差(列间幅度误差)
% % cpp=(channel_error*pi*randn(N,1));cp=cpp*cpp.';%通道间相位误差(列间相位误差)
ch_err_vector = (1+rand(N,1)*channel_error).*exp(2i*pi*rand(N,1)*channel_error); % ####通道误差向量，幅度为以1为均值的高斯分布，相位为0~2π的均匀分布。
if Inn==0 In=ones(N,1);else In=chebwin(N,Inn);end       %根据切比雪夫权值产生数据
if Imm==0 Im=ones(M,1);else Im=chebwin(M,Imm);end       %根据切比雪夫权值产生数据
if Ikk==0 Ik=ones(K,1);else Ik=chebwin(K,Ikk);end       %根据切比雪夫权值产生数据

%====计算不同距离环斜距Rl和高低角phil
Re=6378e3;%地球曲率
rsmax=fun_Rsmax(H);%雷达最大作用斜距
rsu=c/(2*fr);%最大不模糊距离
% 计算（距离模糊）折叠次数L
Nk1 = floor((rsmax-rs)/rsu);                %前向模糊数
Nk2 = floor((rs-H)/rsu);                    %后向模糊数
rsk1 = rs + (0:Nk1)*rsu;                    %前向模糊斜距m
rsk2 = rs - (Nk2:-1:1)*rsu;                 %后向模糊斜距m
if isRu==1
    rsk = [rsk2,rsk1];%Rs;%                 %模糊斜距m
elseif isRu == 0
    rsk = rs;%
end
L = length(rsk);
% L=floor((rsmax-H)/rsu)+1;
% rsk=H+(1:L)*rsu;%####到载机高度距离各次折叠对应的实际距离    0:L的设置使得从H距离处开始采样
% phil=asin(H./rsk);%不同折叠的初始俯仰角

R_percell = 150; % 一个距离门两采样点在空间的距离
N_pcell_max = floor((rsu-H)/R_percell); % 不模糊距离范围内有多少个采样点.注：要保证不要超出距离模糊范围
N_pcell = 30; % ####距离门个数。注：不能大于N_pcell_max
%求取反射率因子G需要的计算
if flag_clutter==1
    F1=4*10^(-7)*10^(0.6*(flag_sea+1));B=pi/2;belt_0=2.44*(flag_sea+1)^1.08/57.29;%不同海情下海杂波系数
elseif  flag_clutter==2
    A=0.00126;B=pi/2;belt_0=0.14;           %沙漠情况下系数
elseif  flag_clutter==3
    A=0.004;B=pi/2;belt_0=0.2;belt_c=1;     %农田情况下系数
elseif  flag_clutter==4
    A=0.0126;B=pi/2;belt_0=0.4;belt_c=1;    %丘陵情况下系数
elseif  flag_clutter==5
    A=0.04;B=1.24;belt_0=0.5;belt_c=1;      %高山情况下系数
end
Pt=180e3;%发射功率
Scale=lambda*sqrt(Pt/(4*pi)^3);%G中的Pt,wavelength,(4*pi)^3
%杂波幅度起伏模型
if flag_kind==3
    clutter_mean=3;clutter_var=7;%K 分布时的U和V(要求是整数)
else
    clutter_mean=0.5;clutter_var=1.2;%log-Normal和weibull杂波均值及方差
end

%% 杂波的频谱分布。方法1为按照风速仿真；方法2为按百分比仿真。不管哪种方法，Vm为0则不产生杂波内部运动，非0值则产生杂波内部运动。

%——方法1:按照风速仿真情况如下（仅考虑风速）
theta_3dB = 0.88*lambda/((N-1)*d);      %主瓣水平向3dB宽度(弧度)
Vm = 0; % ### 如果不产生杂波内部运动，需要将该项设为0，如果需要产生ICM，设任意非0值，因为下面强制对delt_v赋值  %地、海杂波时的风速
% % if flag_clutter==1%公式(2.16)
% %     delta_v1=0.101*Vm;    % 海杂波
% % else
% %     delta_v1=0.0066*Vm;  % 地杂波
% % end
% % delta_v2=V*theta_3dB/(2*sqrt(2*log(2)));%公式(2.17)
% % delt_v=sqrt(delta_v1^2+delta_v2^2);%公式(2.18)
% % % delt_v=sqrt(delta_v1^2);%公式(2.18)
delt_v = 0.5; % 该项为杂波内部运动速度，抛弃上面的设置方式，采取直接强制赋值，杂波速度标准差 0.5m/s
delt_f = 2*delt_v/lambda;%公式(2.15)

%——方法2:按照杂波起伏百分比仿真情况
% % % theta_3dB=0.88*lemda/((N-1)*d);%主瓣水平向3dB宽度(弧度)
% % % delt_f=0.05*2*V/lemda; % 前面的0.05意味着杂波内部运动速度为飞机速度的0.05

f_3db = 2*sqrt(2*log(2))*delt_f;
deltf = 0.6*f_3db;      % 内部运动速度对应的频谱扰动范围
% sk = delt_f*exp(-(deltf*((1:5)-3)).^2/(2*delt_f^2))/sqrt(2*pi*delt_f^2);%公式(2.19)
sk = exp(-(deltf*((1:5)-3)).^2/(2*delt_f^2));%&& 修改了（按J.Ward报告）
Delt = 2*pi*deltf/fr;
sk_sqrt = sqrt(sk);

%% 噪声功率
kk=1.38e-23;%波兹曼常数(J/KHz)
T0=290;%温度
Bn=70e6;%接收机带宽 %% 你确定这样设置没错？
Fn=3.5;%接收机噪声系数dB
Pn=kk*T0*Bn*Fn;%噪声功率
Rorigin=25000;%初始作用距离
f0=c/lambda/10^9;%工作频率GHz
u=sqrt(f0)/4.7;
ops=0;% 运算计数
Sys_DOF = Q*N*K;   % 系统自由度
tic
%% 杂波产生程序
% for ncell=1:N_pcell%第一层循环：不同距离门
for ncell=1:Num_TrainData%第一层循环：不同距离门
    ops=ops+1
    Rcell=Rorigin+(ncell-1)*R_percell;%Rcell为不同距离门对应的不同距离
    clutter_ceho_temp = zeros(Sys_DOF,1);                        % && 用于存储每个距离门的杂波数据
    clutter_covariance_temp = zeros(Sys_DOF);                  % && 用于存储每个距离门的杂波协方差矩阵
    for l=1:L%第二层循环：不同混叠距离
        Rlc=rsk(l)+Rcell; % Rlc为不同距离门不同距离环对应的距离           Rl 里已经加了载机高度 H
        if  Rlc>=H && Rlc<=rsmax
            sin_thetag=H./Rlc-(Rlc.^2-H^2)./(2*Re.*Rlc);%sin(theta_g)
            if sin_thetag>0 && (1-sin_thetag)>(1e-10)
                theta_g=atan(sin_thetag/sqrt(1-sin_thetag^2+(1e-20)));%擦地角
                sin_phil=H/Rlc+(Rlc^2-H^2)/(2*Rlc*Re);%不同距离时的俯仰角的正弦sin(phe)
                phi_l=atan(sin_phil/sqrt(1-sin_phil^2+(1e-20)));%俯仰角
                cos_phil=sqrt(1-sin_phil^2);%俯仰角的余弦
                if flag_clutter==1
                    theta_c=asin(lambda/(4*pi*(0.025+0.046*(flag_sea^1.72))));
                    if theta_g<=theta_c
                        F3=(theta_g/theta_c)^1.9;
                    else
                        F3=1;
                    end
                    delt0=F1*F3*sin(theta_g)/lambda+u*cot(belt_0)^2*exp(-tan(B-theta_g)^2/tan(belt_0)^2);
                elseif flag_clutter==2
                    theta_c=asin(lambda/(4*pi*(9.3*belt_0^2.2)));
                    if theta_g<=theta_c
                        F3=theta_g/theta_c;
                    else
                        F3=1;
                    end
                    delt0=A*F3*sin(theta_g)/lambda+u*cot(belt_0)^2*exp(-tan(B-theta_g)^2/tan(belt_0)^2);
                else
                    delt0=A*sin(theta_g)/lambda+u*cot(belt_0)^2*exp(-tan(B-theta_g)^2/tan(belt_0)^2);
                end
                sr=Rlc*theta_3dB*R_percell/sqrt(2);
                deltc=delt0*sr;%杂波单元雷达截面积
                for ii=1:N_theta%第三层循环：不同方位角
                    if  flag_kind==1
                        al=lognrnd(clutter_mean,clutter_var);
                    elseif  flag_kind==2
                        al=weibell(clutter_var,clutter_mean,0);
                    elseif  flag_kind==3
                        al=-1.9045;%kk_gamm(clutter_mean,clutter_var,0);
                    end
                    if Vm==0
                        X=ones(1,K);
                    else
                        X=cell_f(sk_sqrt,K,Delt);
                    end
                    theta=(ii-1)*(pi/N_theta);
                    ws=2*pi*d*cos(theta)*cos_phil/lambda; % 空域频率
                    wt=4*pi*V*crabM*cos(theta+thetap+crabA)*cos_phil/(lambda*fr); % 时域频率
                    St=exp(1j*wt*(0:K-1)).';
                    Ssr=exp(1j*ws*(0:N-1)).';
                    Sst=exp(1j*ws*gamma*(0:M-1)).';
                    if xinhao==1%%相控阵  gamma=1
                        Us=[1;zeros(M-1,1)];
                        Ut=ones(M,1);
                    elseif xinhao==2%MIMO,gamma>=1
                        Us=eye(M);
                        Ut=eye(M);
                    elseif xinhao==3
                        Us=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
                        Ut=[0.7 0 0 0;0.7 0 0 0;0 0.7 0 0;0 0.7 0 0;0 0 0.7 0;0 0 0.7 0;0 0 0 0.7;0 0 0 0.7];
                    elseif xinhao==4
                        Us=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
                        Ut=[1 0 0 0;1 0 0 0;0 1 0 0;0 1 0 0;0 0 1 0;0 0 1 0;0 0 0 1;0 0 0 1];
                    end
                    %C=Scale*lognrnd(0,1.2,1,1)/Rlc^2*kron(Us.'*Rs*Us*(Ut.')*(Im.*Sst),kron(St,Ssr));
                    C=Scale*sqrt(deltc)/Rlc^2*kron(Us.'*Rs.'*Us*Ut.'*Sst,kron(St,Ssr));
                    C = C/norm(C);
                    clutter_ceho_temp = clutter_ceho_temp+C;                               % && 用于存储每个距离门的杂波数据
                    clutter_covariance_temp = clutter_covariance_temp+C*C';         % && 用于存储每个距离门的杂波协方差矩阵
                end
            end
        end
    end
    clutter_data(:,ncell) = clutter_ceho_temp; % 产生的杂波数据 没加噪声
    clutter_covariance{ncell} = clutter_covariance_temp; % 产生的杂波协方差矩阵
    clear  clutter_ceho_temp   clutter_covariance_temp;
end
toc
parameter=[fr M N K lambda V H el0 az0 crabA thetap];
% save clutterE.mat clutter_data parameter clutter_covariance
ch_num = N; % 阵元数
p_num = K; % 脉冲数
inf_str = '内阵元外脉冲';
% save SL_ULA_8C_8P_300R_IDEAL.mat    clutter_echo    clutter_covariance    inf_str
% save SL_ULA_8C_8P_300R_5.0ICM.mat    clutter_echo    clutter_covariance    inf_str
% save SL_ULA_8C_8P_300R_百分之五ICM.mat    clutter_echo    clutter_covariance    inf_str
% save SL_ULA_8C_8P_300R_0.05CE.mat    ch_err_vector    clutter_echo    clutter_covariance    inf_str

%计算协方差特征值D并绘图
CNR_dB = 60;
CNR = 10^(CNR_dB/10);
R_real = clutter_covariance{1};
%R_real =clutter_echo*clutter_echo'/100;

R = CNR*R_real./sum(eig(R_real)/Sys_DOF);
R = R+1*eye(Sys_DOF);

P.R = clutter_covariance;
P.TrainData = clutter_data;
P.lambda = lambda;
P.Nrow = N;    %每行阵元数
P.Ncol = M;    %每列阵元数
P.Np = K;      %脉冲数      
P.Us = Us;      
P.Ut = Ut; 
P.Rs = Rs;
end

