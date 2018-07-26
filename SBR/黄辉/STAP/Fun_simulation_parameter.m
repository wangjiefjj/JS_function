%function Fun_simulation_parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%环境参数%%%%%%%%%%%%%%%%%%%%%%
entironment_parameter.C=3*10^8;               %光速（m/s)
entironment_parameter.re=8490000;                %地球曲率半径（m)
entironment_parameter.Vm=6;                   %风速（m/s)
entironment_parameter.scene_parameter=5;      %场景参数：1 海杂波 2 沙漠 3 农田 4 丘陵 5 高山 
entironment_parameter.ss=2;                   %海情参数：1 2 3 4 5级海情
%%%%%%%%%%%%%%%%%%%%%%载机参数%%%%%%%%%%%%%%%%%%%%%%%
airborne_parameter.H=6000;                %载机高度(m)
airborne_parameter.V=350;                 %载机速度（m/s)
%%%%%%%%%%%%%%%%%%%%%%信号参数%%%%%%%%%%%%%%%%%%%%%%%
signal_parameter.fz=1.5;           %载频（GHz）
signal_parameter.lambda=0.3/signal_parameter.fz;              %波长(m)
signal_parameter.Bs=0.8*10^6;              %信号带宽
signal_parameter.minu=40*10^(-6);         %脉宽（s）
signal_parameter.fs=1*10^6;              %采样率(次/s)
signal_parameter.fr=4*airborne_parameter.V/signal_parameter.lambda;            %脉冲重复频率（Hz),8*airborne_parameter.V/signal_parameter.lambda
signal_parameter.pulse_num=16;           %脉冲数
signal_parameter.D_zk=signal_parameter.minu*signal_parameter.fr;             %占空比
signal_parameter.D=signal_parameter.Bs*signal_parameter.minu;             %脉压比
signal_parameter.delta_r=entironment_parameter.C/signal_parameter.fs/2;        %距离分辨单元（m）
%signal_parameter.delta_r=150;           %距离分辨单元(km)
%%%%%%%%%%%%%%%%%%%%%阵列参数%%%%%%%%%%%%%%%%%%%%%%%%
natenna_parameter.arry_distance=0.5*signal_parameter.lambda;  %阵元间距
%natenna_parameter.arry_distance=0.5;  %阵元间距
natenna_parameter.arry_col_num=16;                     %阵元列数
natenna_parameter.Wq_azimuth_db=35;          %行加权（db）
natenna_parameter.arry_row_num=1;                     %阵元行数
natenna_parameter.Wq_pitch_db=35;           %列加权（db)

%%%%%%%%%%%%%%%%%%%%%%雷达参数%%%%%%%%%%%%%%%%%%%%%%%
Pt=4.8*10^6;                            %雷达发射的峰值功率(W)
Gt0=1;                                  %雷达发射天线增益
Gr0=1;                                  %雷达接收天线增益
kama=10^(8/10);                                 %系统损耗因子
radar_parameter.r_min=6e3;               %雷达初始作用距离(Km)
radar_parameter.r_max=sqrt(2*entironment_parameter.re*airborne_parameter.H);     %雷达直视距离（m)
radar_parameter.cita_3db=pi/180;             %反射单元角宽(弧度)
%有效接收功率密度:P_r=D_zk*Pt*Gt*Gr*lambda^2/((4*pi)^3*r^4*Kama) r为斜距
radar_parameter.P_r0=Pt*Gt0*Gr0*signal_parameter.lambda^2/((4*pi)^3*kama);
%%%%%%%%%%%%%%%%%%%%%%噪声参数%%%%%%%%%%%%%%%%%%%%%%%
K_borzman=1.38*10^(-23);                %波尔兹曼常数
T_z=290;                                %噪声温度(K)
Bn=0.8*10^6;                             %噪声带宽(Hz)
Fn=4.5;                                 %噪声系数(db)
noise_parameter.Pn=K_borzman*T_z*Bn*10^(Fn/10); %噪声功率(W)
%%%%%%%%%%%%%%%%%%%杂波的幅度起伏模型%%%%%%%%%%%%%%%%%%%%%%
ampfluct_parameter.model=3;              %杂波的幅度起伏分布：1为Log-normal分布，2为Weibull分布 , 3为K分布,K 分布时的U和V(要求是整数) 
ampfluct_parameter.miu=0.2;              %Log-normal分布参数
ampfluct_parameter.sigma=1.2;           %Log-normal分布参数
ampfluct_parameter.b=1.2;                %Weibull分布参数
ampfluct_parameter.c=1.8;                %Weibull分布参数
ampfluct_parameter.u=3;                %K分布参数
ampfluct_parameter.v=7;                %K分布参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








