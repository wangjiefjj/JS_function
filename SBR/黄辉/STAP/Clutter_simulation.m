% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %*******************  主 程 序  *********************%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
isRotation=0;
isRu = 0;
%% 生成仿真参数
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %               杂波仿真模型、参数
% %%%%%%%%%%%%%%%%%%%%%%%环境参数%%%%%%%%%%%%%%%%%%%%%%
entironment_parameter.C=3*10^8;               %光速（m/s)
entironment_parameter.re=8490000;             %地球曲率半径（m)
entironment_parameter.Vm=6;                   %风速（m/s)
entironment_parameter.scene_parameter=5;      %场景参数：1 海杂波 2 沙漠 3 农田 4 丘陵 5 高山 
entironment_parameter.ss=2;                   %海情参数：1 2 3 4 5级海情
%%%%%%%%%%%%%%%%%%%%%%载机参数%%%%%%%%%%%%%%%%%%%%%%%
u=398600.8e9;
SBR_parameter.H=800e3;                %载机高度(m)
re=6378e3;                                 %地球半径
SBR_parameter.V=fun_Vp(SBR_parameter.H);%sqrt(u/(airborne_parameter.H+re))           %平台速度（m/s)
%%%%%%%%%%%%%%%%%%%%%%信号参数%%%%%%%%%%%%%%%%%%%%%%%
signal_parameter.fz=1.25;           %载频（GHz）
signal_parameter.lambda=0.3/signal_parameter.fz;              %波长(m)
signal_parameter.Bs=0.8*10^6;              %信号带宽
signal_parameter.minu=40*10^(-6);         %脉宽（s）
signal_parameter.fs=1*10^6;              %采样率(次/s)
% signal_parameter.fr=4*airborne_parameter.V/signal_parameter.lambda;            %脉冲重复频率（Hz),
signal_parameter.fr=5000;
signal_parameter.pulse_num=16;           %脉冲数
signal_parameter.D_zk=signal_parameter.minu*signal_parameter.fr;             %占空比
signal_parameter.D=signal_parameter.Bs*signal_parameter.minu;             %脉压比
signal_parameter.delta_r=entironment_parameter.C/signal_parameter.fs/2;        %距离分辨单元（m）
%signal_parameter.delta_r=150;           %距离分辨单元(km)
%%%%%%%%%%%%%%%%%%%%%阵列参数%%%%%%%%%%%%%%%%%%%%%%%%
natenna_parameter.arry_distance=0.5*signal_parameter.lambda;  %阵元间距
%natenna_parameter.arry_distance=0.5;  %阵元间距
natenna_parameter.arry_col_num=8;                     %阵元列数
natenna_parameter.Wq_azimuth_db=35;          %行加权（db）
natenna_parameter.arry_row_num=8;                     %阵元行数
natenna_parameter.Wq_pitch_db=35;           %列加权（db)
%%%%%%%%%%%%%%%%%%%%%%雷达参数%%%%%%%%%%%%%%%%%%%%%%%
Pt=4.8*10^8;                            %雷达发射的峰值功率(W)
Gt0=1;                                  %雷达发射天线增益
Gr0=1;                                  %雷达接收天线增益
kama=10^(8/10);                                 %系统损耗因子
radar_parameter.r_min=6e3;               %雷达初始作用距离(Km)
radar_parameter.r_max=sqrt(2*entironment_parameter.re*SBR_parameter.H);     %雷达直视距离（m)
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
ampfluct_parameter.model=2;              %杂波的幅度起伏分布：1为Log-normal分布，2为Weibull分布 
ampfluct_parameter.miu=0.2;              %Log-normal分布参数
ampfluct_parameter.sigma=1.2;           %Log-normal分布参数
ampfluct_parameter.b=1.2;                %Weibull分布参数
ampfluct_parameter.c=1.8;                %Weibull分布参数
ampfluct_parameter.u=3;                %K分布参数
ampfluct_parameter.v=7;                %K分布参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 天基雷达特有参数
SBR_parameter.H=800e3;                %天基高度(m)
Re=6378e3;                                 %地球半径
u=398600.8e9;                              %地心引力常数
alpha1=0/180*pi;%星下点纬度
eta_i=45/180*pi;%轨道倾角
Ve=465.1;%地球赤道速度
Deta=Ve/SBR_parameter.V*(1+SBR_parameter.H/Re);
if isRotation==1
   phi_c=fun_CrabAngle(alpha1/pi*180,eta_i/pi*180,SBR_parameter.H); %单基地偏航角 圆轨道
   rho_c=fun_CrabMagnitude(alpha1/pi*180,eta_i/pi*180,SBR_parameter.H);    %偏航幅度 
else
    phi_c=0;
    rho_c=1;
end

%% 载入参数
re=entironment_parameter.re;
C = entironment_parameter.C;
H=SBR_parameter.H;
V=SBR_parameter.V;
r_max=radar_parameter.r_max;
r_min=radar_parameter.r_min;
P_r0=radar_parameter.P_r0;
delta_r=signal_parameter.delta_r;
pulse_num=signal_parameter.pulse_num;
fr=signal_parameter.fr;
lambda=signal_parameter.lambda;
arry_distance=natenna_parameter.arry_distance;
arry_col_num=natenna_parameter.arry_col_num;
arry_row_num=natenna_parameter.arry_row_num;

%% 设置仿真条件%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

natenna_angle=0*pi/180;       %天线斜侧角
crab_theta=0;                         %偏航角
cita0=pi*(90)/180;                    %发射波束方位角
phi0=pi*(0)/180;                    %发射波束俯仰角
delta_d_cell=0;                         %阵元误差

%%%%%%%%%%%%%%%%%%%杂波的幅度起伏模型参数%%%%%%%%%%%%%%%%%%%%%%
clutter_amp_fluct_flag=0;               %杂波幅度起伏标志位
ampfluct_parameter.model=2;       %杂波的幅度起伏分布：1为Log-normal分布，2为Weibull分布
clutter_freq_fluct_flag=0;              %杂波频率起伏标志位

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_noise =1;
%% 产生杂波%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ru=C/2/fr;                                                          %无模糊距离
r_num=floor(ru/delta_r);                                            %无模糊距离内的距离门数%floor(ru/delta_r)
r_start=1;                                                          %起始距离门floor(r_min/delta_r)
L = ceil((r_max-H)/ru);                                             %距离模糊环数
if isRu==0
    L=1;
end
reflect_cell_num=floor(pi/radar_parameter.cita_3db);                %等距离环划分的反射单元数
cita=pi*(0:reflect_cell_num-1)/reflect_cell_num+natenna_angle;      %杂波反射单元数方位角
Sys_DOF = arry_col_num*arry_row_num*pulse_num;                      %系统自由度
% Clutter=zeros(arry_row_num*arry_col_num*pulse_num,r_num+1);       %杂波存储矩阵(阵元级）
% C_Pt=zeros(1,r_num);                                                                            %杂波峰值功率存储
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for kk=1:1
Clutter=zeros(arry_row_num*arry_col_num*pulse_num,r_num-r_start+1);       %杂波存储矩阵(阵元级）
C_Pt=zeros(1,r_num-r_start+1);                                                                            %杂波峰值功率存储  
R_temp = zeros(arry_row_num*arry_col_num*pulse_num,arry_row_num*arry_col_num*pulse_num);              %杂波协方差矩阵临时存储矩阵
% for n=0:r_num-r_start%n=1:r_num
for n=0:0%n=1:r_num
    rp=(r_start+n)*delta_r;                                                                                                                                        %杂波反射单元径向距离
    
    Clutter_temp=zeros(arry_row_num*arry_col_num,pulse_num);                                                                        %杂波临时存储矩阵(阵元级）
    Clutter_temp_ch=zeros(arry_col_num,pulse_num);                                                                                           %杂波临时存储矩阵(通道级）
    
    
    for l = 0:L-1%L-1
        l
        r=H+rp+ru*l;                                                                                                                                                     %距离模糊环径向距离
        phi=asin(H/r+(r.^2-H^2)/(2*re*r));                                                                                                             %俯仰角
        
%         F_t=function_Ft(cita,phi,cita0,phi0);                                                                                                                %发射方向图
        [F_t(l+1,:),F_r]=Fun_trans_pattern(cita0,phi0,natenna_angle,crab_theta,cita,phi,lambda,natenna_parameter);                                  %发射和接收方向图
        P_r=P_r0/r^4;                                                                                                                                              %反射功率
        Sigma=Fun_clutter_reflect_model(r,phi,SBR_parameter,radar_parameter,signal_parameter,entironment_parameter);   %雷达等效截面积
        clutter_amp=sqrt(P_r*Sigma).*F_t(l+1,:);                                                                                                                                                %等距离环杂波信号的幅度
        
        clutter_amp_fluct=Fun_clutter_amp_fluct(clutter_amp_fluct_flag,reflect_cell_num,ampfluct_parameter,signal_parameter);   %等距离环各单元格杂波信号的幅度起伏
        clutter_freq_fluct=Fun_clutter_freq_fluct(clutter_freq_fluct_flag,reflect_cell_num,entironment_parameter,SBR_parameter,signal_parameter,radar_parameter);    %杂波信号的频率起伏
        
        Phase_rand=exp(1j*pi*(1-2*rand(1,reflect_cell_num)));                                                                                %等距离环各单元格杂波信号的随机相位
        
        %F_t=ones(1,reflect_cell_num);
        %arry_cell_err=(1+delta_d_cell*randn(arry_row_num,arry_col_num)).*exp(j*pi*delta_d_cell*(unifrnd(-pi,pi,arry_row_num,arry_col_num)));%产生阵元误差
        
        C_0=(clutter_amp.*clutter_amp_fluct.*Phase_rand).'*ones(1,pulse_num).*clutter_freq_fluct;                                     %基准阵元接收到pulse_num脉冲的杂波数据（无时频）
        t_frq=2*rho_c*cos(phi)*cos(cita+phi_c)*V/lambda;         %时频
        t_ph=exp(1j*2*pi*(t_frq.'*(0:pulse_num-1)/fr));                                                                                            %时频相位
        C_z=C_0.*t_ph;                                                                                                                                          %基准阵元接收到pulse_num脉冲的杂波数据（含时频）
        C_Pt(n+1)=clutter_amp*clutter_amp'+C_Pt(n+1);
        
        %% 无通道合成数据
        for m=1:reflect_cell_num
            s_frq_azimuth=cos(phi)*cos(cita(m)-natenna_angle)*arry_distance/lambda;                          %方位维空频
            s_frq_pitch=sin(phi+crab_theta)*arry_distance/lambda;                                %俯仰维空频
            s_ph_azimuth=exp(-1j*2*pi*s_frq_azimuth.'*(0:arry_col_num-1));                %空频方位维相位
            s_ph_pitch=exp(-1j*2*pi*s_frq_pitch.'*(0:arry_row_num-1));                        %空频俯仰维相位
            arry_cell_err=(1+delta_d_cell*randn(arry_row_num*arry_col_num,1)).*exp(j*pi*delta_d_cell*(unifrnd(-1,1,arry_row_num*arry_col_num,1)));%产生阵元误差
            s_ph(m,:)=kron(s_ph_pitch,s_ph_azimuth).*arry_cell_err.';                                                      %合成的空频相位
            %%计算真实的协方差矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if (r_start+n)==ceil((r_num+r_start)/2)
            X_temp = kron((C_0(m,:).*t_ph(m,:)).',s_ph(m,:).');
            R_temp = X_temp*X_temp' + R_temp;
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        Clutterz=s_ph.'*(C_z);                                                                                  %所有阵元接收到当前距离环的pulse_num脉冲的杂波数据
        Clutter_temp=Clutter_temp+Clutterz;                                                       %对所有距离模糊环杂波数据进行累积
        
    end
    
    %% 产生噪声%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flag_noise == 1
        Z =randn(Sys_DOF,1);
        gauss_noise = sqrt(noise_parameter.Pn)*Z.*exp(1j*unifrnd(-pi,pi,Sys_DOF,1));%产生高斯噪声
    else
        gauss_noise = zeros(Sys_DOF,1);
    end
    
    %% 添加噪声    
    Clutter(:,n+1)= reshape(Clutter_temp,arry_row_num*arry_col_num*pulse_num,1)+gauss_noise;    
    Clutter_covariance(:,:,n+1)=R_temp+noise_parameter.Pn*eye(Sys_DOF,Sys_DOF);
end
clutter_data{kk}.clutter= Clutter;
clutter_data{kk}.CNR_real=10*log10(C_Pt/noise_parameter.Pn); 
clutter_data{kk}.clutter_covariance=Clutter_covariance;%R_temp+noise_parameter.Pn*eye(Sys_DOF,Sys_DOF);%;
end
Pn=noise_parameter.Pn;
parameter=[fr arry_row_num arry_col_num pulse_num lambda V arry_distance phi0 cita0 Pn];
% save clutter_data16x16.mat clutter_data parameter
toc

Sys_DOF=arry_col_num*arry_row_num*pulse_num;
% clutter_data= Clutter;


% CNR_real=10*log10(C_Pt/noise_parameter.Pn);
% noise_covariance=noise_parameter.Pn*eye(Sys_DOF,Sys_DOF);
% save clutter_data.mat clutter_data parameter CNR_real
% save clutter_covariance.mat clutter_covariance noise_covariance


% %% 杂波功率谱分析
% 
% % 计算协方差特征值D并绘图
R =clutter_data{1}.clutter_covariance(:,:,1);
% R =Clutter*Clutter'/150;
CNR_dB = 60;
CNR = 10^(CNR_dB/10);
%R_real =Clutter_ch*Clutter_ch'/30;
%R_real =Clutter*Clutter'/30;
R_real = R;
R = CNR*R_real./sum(eig(R_real)/Sys_DOF);
R = R+1*eye(Sys_DOF);

DD= eig(R);
DD = sort(abs(DD),'descend');
D_db = db(DD)/2;
figure;plot(D_db','*-');grid on;title('杂波协方差矩阵的特征谱');%-max(D_db)
% % 
% % % % 计算杂波功率谱P并绘图
% % Nfs=180;
% % fs = (-Nfs:1:Nfs)/Nfs;
% % Nfd=180;
% % fd = (-Nfd:1:Nfd)/Nfd;
% % fd=2*cos(phi0)*cos(pi*fs+natenna_angle+crab_theta)*V/lambda/fr;
% % R_inv=inv(R);
% % for p=1:length(fs)
% %     s_ph_azimuth=exp(-j*pi*fs(p)*(0:arry_col_num-1));                         %空频方位维相位
% %     s_ph_pitch=exp(-j*pi*fs(p)*(0:arry_row_num-1));                                 %空频俯仰维相位
% %     Ss=kron(s_ph_pitch,s_ph_azimuth);                                                                   %合成的空频相位
% %     for q=1:length(fd)
% %         St=exp(j*pi*(0:pulse_num-1).'*fd(q));
% %         S=kron(St,Ss.');
% %         S = S/norm(S);
% %         P(p,q)=1/abs(S'*R_inv*S);
% %     end
% % end
% % P_db = db(P)/2;
% % P_db =P_db -max(max(P_db ));
% % P_db =P_db ;
% % figure; %surf(P_db);shading interp;
% % axis([-0.5,0.5,-0.5,0.5,-60,0]);
% % angle=180*acos(fs)/pi;
% %  mesh(fd,fs,P_db);
% % surf(fd,fs,P_db);shading interp;lighting gouraud;colorbar;axis tight;
% % title('杂波功率谱');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB');  %绘制杂波功率谱图mesh(fd,fs,P_db);cos(\psi)
% % figure;contour(fd,fs,P_db,20);   %绘制杂波功率谱等高线图
% % title('杂波功率谱等高线');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB');  
% %
% %

