clc;clear;close all;
% load clutter_data16x16.mat
% load clutter_data_target.mat
%load target.mat;
% Cr1 = clutter_data{1}.clutter;
% load clutter_data.mat
% Cr2 = clutter_data{1}.clutter;
% load  Channel_12_Pulse_16_ClutterData_Mountain.mat
% load  Channel_12_Pulse_16_ClutterData_Sea.mat
% parameter=[fr arry_row_num arry_col_num pulse_num lambda V H phi0 cita0 crab_theta natenna_angle];
%% 系统参数
V = parameter(6);%飞机速度
N = parameter(3);%每行阵元数
M = parameter(2);%每列阵元数
K = parameter(4);%脉冲数
lemda = parameter(5);%波长
d = parameter(7);%阵元间距
fr = parameter(1);%脉冲重复频率
phi0 = parameter(8);%波束指向俯仰角
theta0 = parameter(9);%波束指向方位角
% thetap = parameter(11);%阵面和航迹之间夹角/sqrt(N*K)
%  Cr = clutter_data{1}.clutter;%杂波数据矢量，摆放顺序为俯仰、方位和脉冲。
 Cov=clutter_data{1}.clutter_covariance(:,:,1);%杂波数据协方差矩阵
 L=3000;
Cr=Cov^(0.5)*((randn(N*K,L)).*exp(j*unifrnd(-pi,pi,N*K,L)));
ERR=inv(Cov)*(Cr*Cr'/L);
 %Cr = target;%杂波数据矢量，摆放顺序为俯仰、方位和脉冲。
%Cov=R_cov;%杂波数据协方差矩阵
CNR_dB =60;%杂噪比dB
CNR = 10^(CNR_dB/10);%杂噪比
 channel_error = 0;%通道误差
%% 空域通道合成
ws_azimuth_sig = 2*pi*d*cos(phi0)*cos(theta0)/lemda;
ws_pitch_sig = 2*pi*d*sin(phi0)/lemda;
Ss_azimuth_sig = exp(j*(0:N-1).'*ws_azimuth_sig);%目标向方位空域导向矢量
Ss_pitch_sig = exp(j*(0:M-1).'*ws_pitch_sig);%目标向俯仰空域导向矢量
Ss_sig = kron(Ss_pitch_sig,Ss_azimuth_sig);%目标向空域导向矢量
 
Ns = N;%方位维合成子阵个数
Ms = M;%俯仰维合成子阵个数

%% 子阵合成
T_synthesize= function_sub_arry_synthesize(N,Ns,M,Ms,K,Ss_pitch_sig,Ss_azimuth_sig,channel_error);
Cr_T =T_synthesize*Cr;
Ss_azimuth_syn = exp(j*(0:Ns-1).'*ws_azimuth_sig);%目标向方位空域导向矢量
Ss_pitch_syn = exp(j*(0:Ms-1).'*ws_pitch_sig);%目标向俯仰空域导向矢量
Ss_sig_syn= kron(Ss_pitch_syn,Ss_azimuth_syn);%目标向空域导向矢量
Sys_DOF=Ns*Ms*K;
%% STAP处理
% [theta_a,theta_e,fd_s]=function_partion_APF(Ns,Ms,K,d,V,lemda,fr);
fd_s=linspace(-1,1,2*K+1);
Data=Cr_T;      %Data：当前景数据
fd= fd_s;     %fd:需检测得多普勒通道
% Rangegate_now=500;      %Rangegate_now:当前待检测距离门
Rangegate_end=length(Data);      %Rangegate_end:最大距离门length(Data)
Ss=Ss_sig_syn;      %Ss:待检测单元空域导向矢量
Wq=chebwin(K,40);        %Wq:时域通道加权 
Ws=chebwin(Ns*Ms,40);        %Ws:空域通道加权
 %% OPT
N_TrainData=2*Ns*Ms*K;      %N_TrainData：训练样本数
  for Rangegate_now=25;
 %for     Rangegate_now=500;
%  R_real=T_synthesize*clutter_covariance{Rangegate_now}*T_synthesize';    %R_real:待检测单元真实的R阵
%  R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
 R_real=T_synthesize*Cov*T_synthesize'; 
STAP_OPT= function_STAP_OPT(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,R_real,Ss,K,CNR);
 STAP_OPT_Output(:,Rangegate_now)=STAP_OPT.Output;
 end

%% JDL
N_timechannel=3;%N_timechannel:选取时域通道数
N_spacechannel=3;%N_spacechannel:选取空域通道数
space_F=[-floor(N_spacechannel/2):1:floor(N_spacechannel/2)]/Ns;
ws_azimuth_sig = 2*pi*d*(cos(theta0)*cos(phi0)+space_F)/lemda;
ws_pitch_sig = 2*pi*d*sin(phi0)/lemda;
Ss_azimuth_sig = exp(j*(0:Ns-1).'*ws_azimuth_sig);%目标向方位空域导向矢量
Ss_pitch_sig = exp(j*(0:Ms-1).'*ws_pitch_sig);%目标向俯仰空域导向矢量
Ts = diag(Ws)*kron(Ss_azimuth_sig,Ss_pitch_sig);%Ts:空域通道导向矢量矩阵
%Wq=ones(K,1);
N_TrainData=3*N_timechannel*N_spacechannel;

% for Rangegate_now=2:Rangegate_end
for     Rangegate_now=25;    
% R_real=T_synthesize*Cov(:,:,Rangegate_now)*T_synthesize';    %R_real:待检测单元真实的R阵
R_real=T_synthesize*Cov*T_synthesize'; 
% R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
STAP_JDL= function_STAP_JDL(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,R_real,Ts,Ss,Wq,K,CNR);
STAP_JDL_Output(:,Rangegate_now)=STAP_JDL.Output;
end
%% 3DT

arry_num=Ns*Ms;
N_timechannel=3;%N_timechannel:选取时域通道数
N_TrainData=2*N_timechannel*arry_num;

for Rangegate_now=25
%for     Rangegate_now=1;
% R_real=T_synthesize*Cov(:,:,Rangegate_now)*T_synthesize';    %R_real:待检测单元真实的R阵  
R_real=T_synthesize*Cov*T_synthesize'; 
% R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
STAP_MCAP= function_STAP_MCAP(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,R_real,arry_num,Ss,Wq,K,CNR);
STAP_MCAP_Output(:,Rangegate_now)=STAP_MCAP.Output;
end

% %% GMB
% N_timechannel=7;%N_timechannel:选取时域通道数
% N_spacechannel=9;%N_spacechannel:选取空域通道数
% space_F=[-floor(N_spacechannel/2):1:floor(N_spacechannel/2)]/Ns;
% ws_azimuth_sig = 2*pi*d*(cos(theta0)*cos(phi0)+space_F)/lemda;
% ws_pitch_sig = 2*pi*d*sin(phi0)/lemda;
% Ss_azimuth_sig = exp(j*(0:Ns-1).'*ws_azimuth_sig);%目标向方位空域导向矢量
% Ss_pitch_sig = exp(j*(0:Ms-1).'*ws_pitch_sig);%目标向俯仰空域导向矢量
% Ts =diag(Ws)* kron(Ss_azimuth_sig,Ss_pitch_sig);%Ts:空域通道导向矢量矩阵
% %Wq=ones(K,1);
% N_TrainData=10*(N_timechannel+N_spacechannel-1);
% 
% %for Rangegate_now=2:Rangegate_end
% for     Rangegate_now=25;
% % R_real=T_synthesize*Cov(:,:,Rangegate_now)*T_synthesize';    %R_real:待检测单元真实的R阵
% R_real=T_synthesize*Cov*T_synthesize'; 
% R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
% 
% STAP_GMB= function_STAP_GMB(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,N_spacechannel,R_real,Ts,Ss,Wq,K,CNR);
% STAP_GMB_Output(:,Rangegate_now)=STAP_GMB.Output;
% end

%%  F$A
%%  A$F



%% 画改善因子
figure(1)
plot( fd,0.5*db(abs(STAP_OPT.IF_OPT)),'m--o','markersize',15,'LineWidth',4);
hold on
plot( fd,0.5*db(abs(STAP_JDL.IF_JDL)),'r--d','markersize',15,'LineWidth',4);
hold on
plot( fd,0.5*db(abs(STAP_MCAP.IF_MCAP)),'b--x','markersize',15,'LineWidth',4);

% hold on
% plot( 0.5*db(abs(STAP_GMB.IF_GMB)),'-*k');
legend('OPT','JDL','3DT','FontSize',80,'location','SouthEast');
xlabel('2fd/fr','FontSize',40);ylabel('IF/dB','FontSize',40);
set(gca,'FontSize',40);
grid on
%% 用于画二维响应
% Nfs=180;
% Nfd=ceil(Nfs*fr/(4*V/lemda));
% fd = linspace(-1,1,Nfd);
% fs = linspace(-1,1,Nfs);
% theta_a=fs;fd = linspace(-1,1,65);
theta_a=linspace(-1,1,360);
ws_azimuth_response = 2*pi*d*(cos(theta0)*cos(phi0)+theta_a)/lemda;
ws_pitch_response = 2*pi*d*sin(phi0)/lemda;
Ss_azimuth_response = exp(j*(0:Ns-1).'*ws_azimuth_response);%目标向方位空域导向矢量
Ss_pitch_response = exp(j*(0:Ms-1).'*ws_pitch_response);%目标向俯仰空域导向矢量
Ss_sig_response= kron(Ss_azimuth_response,Ss_pitch_response);%目标向空域导向矢量
Nfd = length(fd);
% %%%%%%%%
W_OPT=STAP_OPT.W_OPT(:,6);
for n = 1:Nfd     
    St = exp(j*[0:K-1]'*pi*fd(n));
    S1 = kron(St,Ss_sig_response);
    Respone_OPT(:,n)=db(abs(W_OPT'*S1));
end
figure(2);surf(Respone_OPT);shading interp;lighting gouraud;colorbar;axis tight;
W_JDL=STAP_JDL.W_JDL(:,6);
for n = 1:Nfd     
    St = exp(j*[0:K-1]'*pi*fd(n));
    S2 = kron(St,Ss_sig_response);
    Respone_JDL(:,n)=db(abs(W_JDL'*S2));    
end
figure(3);surf(Respone_JDL);shading interp;lighting gouraud;colorbar;axis tight;
W_MCAP=STAP_MCAP.W_MCAP(:,6);
for n = 1:Nfd     
    St = exp(j*[0:K-1]'*pi*fd(n));
    S3 = kron(St,Ss_sig_response);
    Respone_MCAP(:,n)=db(abs(W_MCAP'*S3));
end
figure(4);surf(Respone_MCAP);shading interp;lighting gouraud;colorbar;axis tight;
% W_GMB=STAP_GMB.W_GMB(:,6);
% for n = 1:Nfd     
%     St = exp(j*[0:K-1]'*pi*fd(n));
%     S4 = kron(St,Ss_sig_response);
%     Respone_GMB(:,n)=db(abs(W_GMB'*S4)); 
% end
% Respone_GMB(find(Respone_GMB<=-110))=-110;
% figure(5);surf(Respone_GMB);shading interp;lighting gouraud;colorbar;axis tight;
%% 杂波抑制前后的距离多普勒图
% %杂波抑制前
% figure(6)
% fd_s=linspace(-1,1,100);
% Data=Cr_T;      %Data：当前景数据
% fd= fd_s;     %fd:待检测多普勒通道
% St = exp(j*[0:K-1]'*pi*fd);%待检测多普勒通道时域导向矢量
% Ss=Ss_sig_syn;%待检测空域导向矢量
% S=kron(St,Ss);
% Output_B=abs(S'*Data);
% surf(db(Output_B));shading interp;lighting gouraud;colorbar;axis tight;

%杂波抑制后
% figure(7)
% surf(db(STAP_OPT_Output));shading interp;lighting gouraud;colorbar;axis tight;

% figure(8);
% surf(db(STAP_JDL_Output));shading interp;lighting gouraud;colorbar;axis tight;
% AAA=STAP_MCAP_Output(:,201);
% AAA(find(AAA<10^(-9)))=10^(-9);
% figure(9);
% plot(db(AAA));
% surf(db(STAP_MCAP_Output));shading interp;lighting gouraud;colorbar;axis tight;

% figure(10);
% surf(db(STAP_GMB_Output));shading interp;lighting gouraud;colorbar;axis tight;
% % 
% 
% 
% 
% %% 画杂波谱
% %计算协方差特征值D并绘图
% Sys_DOF=N*K;
% CNR_dB = 60;
% CNR = 10^(CNR_dB/10);
% %R_real =Clutter_ch*Clutter_ch'/30;
% %R_real =Clutter*Clutter'/30;
% R_real =clutter_covariance{Rangegate_now};
% R = CNR*R_real./sum(eig(R_real)/Sys_DOF);
% R = R+1*eye(Sys_DOF);
% [VV,DD]= eig(R);
% DD = sort(abs(DD),'descend');
% D_db = db(DD)/2;
% figure(11);plot(D_db','*-');grid on;title('杂波协方差矩阵的特征谱');
% 
% %计算杂波功率谱P并绘图
% Nfs=180;
% Nfd=ceil(Nfs*fr/(4*V/lemda));
% fd = linspace(-1,1,Nfd);
% fs = linspace(-1,1,Nfs);
% % receive_frq_pitch=sin(phi0)*arry_distance/lemda;
% % receive_ph_pitch=exp(j*2*pi*receive_frq_pitch.'*(0:arry_row_num-1));
% % % ch_combine_matric=kron(kron(receive_ph_pitch,eye(arry_col_num),ones(arry_col_num));
% % % R_ch=conj(ch_combine_matric)*R*ch_combine_matric.';
% 
% R_inv=inv(R);
% for p=1:length(fs)
%     receive_ph_azimuth=exp(j*pi*(0:N-1).'*fs(p));
%     %     Ss=kron(receive_ph_pitch.',receive_ph_azimuth);
%     Ss=receive_ph_azimuth;
%     for q=1:length(fd)
%         St=exp(j*pi*(0:K-1).'*fd(q));
%         S=kron(St,Ss);
%         S = S/norm(S);
%         P(p,q)=1/abs(S'*R_inv*S);
%     end
% end
% P_db = db(P)/2;
% P_db =P_db -max(max(P_db ));
% %P_db =P_db ;
% figure(12);%surf(P_db);shading interp;
% axis([-0.5,0.5,-0.5,0.5,-60,0]);mesh(fd,fs,P_db);
% title('杂波功率谱');xlabel('2fd/fr');ylabel('cos\psi');zlabel('P/dB');  %绘制杂波功率谱图mesh(fd,fs,P_db);
% figure(13);contour(P_db,20);   %绘制杂波功率谱等高线图
