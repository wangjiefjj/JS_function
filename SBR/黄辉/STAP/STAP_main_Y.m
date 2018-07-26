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
%% ϵͳ����
V = parameter(6);%�ɻ��ٶ�
N = parameter(3);%ÿ����Ԫ��
M = parameter(2);%ÿ����Ԫ��
K = parameter(4);%������
lemda = parameter(5);%����
d = parameter(7);%��Ԫ���
fr = parameter(1);%�����ظ�Ƶ��
phi0 = parameter(8);%����ָ������
theta0 = parameter(9);%����ָ��λ��
% thetap = parameter(11);%����ͺ���֮��н�/sqrt(N*K)
%  Cr = clutter_data{1}.clutter;%�Ӳ�����ʸ�����ڷ�˳��Ϊ��������λ�����塣
 Cov=clutter_data{1}.clutter_covariance(:,:,1);%�Ӳ�����Э�������
 L=3000;
Cr=Cov^(0.5)*((randn(N*K,L)).*exp(j*unifrnd(-pi,pi,N*K,L)));
ERR=inv(Cov)*(Cr*Cr'/L);
 %Cr = target;%�Ӳ�����ʸ�����ڷ�˳��Ϊ��������λ�����塣
%Cov=R_cov;%�Ӳ�����Э�������
CNR_dB =60;%�����dB
CNR = 10^(CNR_dB/10);%�����
 channel_error = 0;%ͨ�����
%% ����ͨ���ϳ�
ws_azimuth_sig = 2*pi*d*cos(phi0)*cos(theta0)/lemda;
ws_pitch_sig = 2*pi*d*sin(phi0)/lemda;
Ss_azimuth_sig = exp(j*(0:N-1).'*ws_azimuth_sig);%Ŀ����λ������ʸ��
Ss_pitch_sig = exp(j*(0:M-1).'*ws_pitch_sig);%Ŀ������������ʸ��
Ss_sig = kron(Ss_pitch_sig,Ss_azimuth_sig);%Ŀ���������ʸ��
 
Ns = N;%��λά�ϳ��������
Ms = M;%����ά�ϳ��������

%% ����ϳ�
T_synthesize= function_sub_arry_synthesize(N,Ns,M,Ms,K,Ss_pitch_sig,Ss_azimuth_sig,channel_error);
Cr_T =T_synthesize*Cr;
Ss_azimuth_syn = exp(j*(0:Ns-1).'*ws_azimuth_sig);%Ŀ����λ������ʸ��
Ss_pitch_syn = exp(j*(0:Ms-1).'*ws_pitch_sig);%Ŀ������������ʸ��
Ss_sig_syn= kron(Ss_pitch_syn,Ss_azimuth_syn);%Ŀ���������ʸ��
Sys_DOF=Ns*Ms*K;
%% STAP����
% [theta_a,theta_e,fd_s]=function_partion_APF(Ns,Ms,K,d,V,lemda,fr);
fd_s=linspace(-1,1,2*K+1);
Data=Cr_T;      %Data����ǰ������
fd= fd_s;     %fd:����ö�����ͨ��
% Rangegate_now=500;      %Rangegate_now:��ǰ����������
Rangegate_end=length(Data);      %Rangegate_end:��������length(Data)
Ss=Ss_sig_syn;      %Ss:����ⵥԪ������ʸ��
Wq=chebwin(K,40);        %Wq:ʱ��ͨ����Ȩ 
Ws=chebwin(Ns*Ms,40);        %Ws:����ͨ����Ȩ
 %% OPT
N_TrainData=2*Ns*Ms*K;      %N_TrainData��ѵ��������
  for Rangegate_now=25;
 %for     Rangegate_now=500;
%  R_real=T_synthesize*clutter_covariance{Rangegate_now}*T_synthesize';    %R_real:����ⵥԪ��ʵ��R��
%  R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
 R_real=T_synthesize*Cov*T_synthesize'; 
STAP_OPT= function_STAP_OPT(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,R_real,Ss,K,CNR);
 STAP_OPT_Output(:,Rangegate_now)=STAP_OPT.Output;
 end

%% JDL
N_timechannel=3;%N_timechannel:ѡȡʱ��ͨ����
N_spacechannel=3;%N_spacechannel:ѡȡ����ͨ����
space_F=[-floor(N_spacechannel/2):1:floor(N_spacechannel/2)]/Ns;
ws_azimuth_sig = 2*pi*d*(cos(theta0)*cos(phi0)+space_F)/lemda;
ws_pitch_sig = 2*pi*d*sin(phi0)/lemda;
Ss_azimuth_sig = exp(j*(0:Ns-1).'*ws_azimuth_sig);%Ŀ����λ������ʸ��
Ss_pitch_sig = exp(j*(0:Ms-1).'*ws_pitch_sig);%Ŀ������������ʸ��
Ts = diag(Ws)*kron(Ss_azimuth_sig,Ss_pitch_sig);%Ts:����ͨ������ʸ������
%Wq=ones(K,1);
N_TrainData=3*N_timechannel*N_spacechannel;

% for Rangegate_now=2:Rangegate_end
for     Rangegate_now=25;    
% R_real=T_synthesize*Cov(:,:,Rangegate_now)*T_synthesize';    %R_real:����ⵥԪ��ʵ��R��
R_real=T_synthesize*Cov*T_synthesize'; 
% R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
STAP_JDL= function_STAP_JDL(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,R_real,Ts,Ss,Wq,K,CNR);
STAP_JDL_Output(:,Rangegate_now)=STAP_JDL.Output;
end
%% 3DT

arry_num=Ns*Ms;
N_timechannel=3;%N_timechannel:ѡȡʱ��ͨ����
N_TrainData=2*N_timechannel*arry_num;

for Rangegate_now=25
%for     Rangegate_now=1;
% R_real=T_synthesize*Cov(:,:,Rangegate_now)*T_synthesize';    %R_real:����ⵥԪ��ʵ��R��  
R_real=T_synthesize*Cov*T_synthesize'; 
% R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
STAP_MCAP= function_STAP_MCAP(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,R_real,arry_num,Ss,Wq,K,CNR);
STAP_MCAP_Output(:,Rangegate_now)=STAP_MCAP.Output;
end

% %% GMB
% N_timechannel=7;%N_timechannel:ѡȡʱ��ͨ����
% N_spacechannel=9;%N_spacechannel:ѡȡ����ͨ����
% space_F=[-floor(N_spacechannel/2):1:floor(N_spacechannel/2)]/Ns;
% ws_azimuth_sig = 2*pi*d*(cos(theta0)*cos(phi0)+space_F)/lemda;
% ws_pitch_sig = 2*pi*d*sin(phi0)/lemda;
% Ss_azimuth_sig = exp(j*(0:Ns-1).'*ws_azimuth_sig);%Ŀ����λ������ʸ��
% Ss_pitch_sig = exp(j*(0:Ms-1).'*ws_pitch_sig);%Ŀ������������ʸ��
% Ts =diag(Ws)* kron(Ss_azimuth_sig,Ss_pitch_sig);%Ts:����ͨ������ʸ������
% %Wq=ones(K,1);
% N_TrainData=10*(N_timechannel+N_spacechannel-1);
% 
% %for Rangegate_now=2:Rangegate_end
% for     Rangegate_now=25;
% % R_real=T_synthesize*Cov(:,:,Rangegate_now)*T_synthesize';    %R_real:����ⵥԪ��ʵ��R��
% R_real=T_synthesize*Cov*T_synthesize'; 
% R_real = CNR*R_real./sum(eig(R_real)/Sys_DOF);R_real = R_real+1*eye(Sys_DOF);
% 
% STAP_GMB= function_STAP_GMB(Data,fd,Rangegate_now,Rangegate_end,N_TrainData,N_timechannel,N_spacechannel,R_real,Ts,Ss,Wq,K,CNR);
% STAP_GMB_Output(:,Rangegate_now)=STAP_GMB.Output;
% end

%%  F$A
%%  A$F



%% ����������
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
%% ���ڻ���ά��Ӧ
% Nfs=180;
% Nfd=ceil(Nfs*fr/(4*V/lemda));
% fd = linspace(-1,1,Nfd);
% fs = linspace(-1,1,Nfs);
% theta_a=fs;fd = linspace(-1,1,65);
theta_a=linspace(-1,1,360);
ws_azimuth_response = 2*pi*d*(cos(theta0)*cos(phi0)+theta_a)/lemda;
ws_pitch_response = 2*pi*d*sin(phi0)/lemda;
Ss_azimuth_response = exp(j*(0:Ns-1).'*ws_azimuth_response);%Ŀ����λ������ʸ��
Ss_pitch_response = exp(j*(0:Ms-1).'*ws_pitch_response);%Ŀ������������ʸ��
Ss_sig_response= kron(Ss_azimuth_response,Ss_pitch_response);%Ŀ���������ʸ��
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
%% �Ӳ�����ǰ��ľ��������ͼ
% %�Ӳ�����ǰ
% figure(6)
% fd_s=linspace(-1,1,100);
% Data=Cr_T;      %Data����ǰ������
% fd= fd_s;     %fd:����������ͨ��
% St = exp(j*[0:K-1]'*pi*fd);%����������ͨ��ʱ����ʸ��
% Ss=Ss_sig_syn;%����������ʸ��
% S=kron(St,Ss);
% Output_B=abs(S'*Data);
% surf(db(Output_B));shading interp;lighting gouraud;colorbar;axis tight;

%�Ӳ����ƺ�
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
% %% ���Ӳ���
% %����Э��������ֵD����ͼ
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
% figure(11);plot(D_db','*-');grid on;title('�Ӳ�Э��������������');
% 
% %�����Ӳ�������P����ͼ
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
% title('�Ӳ�������');xlabel('2fd/fr');ylabel('cos\psi');zlabel('P/dB');  %�����Ӳ�������ͼmesh(fd,fs,P_db);
% figure(13);contour(P_db,20);   %�����Ӳ������׵ȸ���ͼ
