% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %*******************  �� �� ��  *********************%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
isRotation=0;
isRu = 0;
%% ���ɷ������
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %               �Ӳ�����ģ�͡�����
% %%%%%%%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%
entironment_parameter.C=3*10^8;               %���٣�m/s)
entironment_parameter.re=8490000;             %�������ʰ뾶��m)
entironment_parameter.Vm=6;                   %���٣�m/s)
entironment_parameter.scene_parameter=5;      %����������1 ���Ӳ� 2 ɳĮ 3 ũ�� 4 ���� 5 ��ɽ 
entironment_parameter.ss=2;                   %���������1 2 3 4 5������
%%%%%%%%%%%%%%%%%%%%%%�ػ�����%%%%%%%%%%%%%%%%%%%%%%%
u=398600.8e9;
SBR_parameter.H=800e3;                %�ػ��߶�(m)
re=6378e3;                                 %����뾶
SBR_parameter.V=fun_Vp(SBR_parameter.H);%sqrt(u/(airborne_parameter.H+re))           %ƽ̨�ٶȣ�m/s)
%%%%%%%%%%%%%%%%%%%%%%�źŲ���%%%%%%%%%%%%%%%%%%%%%%%
signal_parameter.fz=1.25;           %��Ƶ��GHz��
signal_parameter.lambda=0.3/signal_parameter.fz;              %����(m)
signal_parameter.Bs=0.8*10^6;              %�źŴ���
signal_parameter.minu=40*10^(-6);         %����s��
signal_parameter.fs=1*10^6;              %������(��/s)
% signal_parameter.fr=4*airborne_parameter.V/signal_parameter.lambda;            %�����ظ�Ƶ�ʣ�Hz),
signal_parameter.fr=5000;
signal_parameter.pulse_num=16;           %������
signal_parameter.D_zk=signal_parameter.minu*signal_parameter.fr;             %ռ�ձ�
signal_parameter.D=signal_parameter.Bs*signal_parameter.minu;             %��ѹ��
signal_parameter.delta_r=entironment_parameter.C/signal_parameter.fs/2;        %����ֱ浥Ԫ��m��
%signal_parameter.delta_r=150;           %����ֱ浥Ԫ(km)
%%%%%%%%%%%%%%%%%%%%%���в���%%%%%%%%%%%%%%%%%%%%%%%%
natenna_parameter.arry_distance=0.5*signal_parameter.lambda;  %��Ԫ���
%natenna_parameter.arry_distance=0.5;  %��Ԫ���
natenna_parameter.arry_col_num=8;                     %��Ԫ����
natenna_parameter.Wq_azimuth_db=35;          %�м�Ȩ��db��
natenna_parameter.arry_row_num=8;                     %��Ԫ����
natenna_parameter.Wq_pitch_db=35;           %�м�Ȩ��db)
%%%%%%%%%%%%%%%%%%%%%%�״����%%%%%%%%%%%%%%%%%%%%%%%
Pt=4.8*10^8;                            %�״﷢��ķ�ֵ����(W)
Gt0=1;                                  %�״﷢����������
Gr0=1;                                  %�״������������
kama=10^(8/10);                                 %ϵͳ�������
radar_parameter.r_min=6e3;               %�״��ʼ���þ���(Km)
radar_parameter.r_max=sqrt(2*entironment_parameter.re*SBR_parameter.H);     %�״�ֱ�Ӿ��루m)
radar_parameter.cita_3db=pi/180;             %���䵥Ԫ�ǿ�(����)
%��Ч���չ����ܶ�:P_r=D_zk*Pt*Gt*Gr*lambda^2/((4*pi)^3*r^4*Kama) rΪб��
radar_parameter.P_r0=Pt*Gt0*Gr0*signal_parameter.lambda^2/((4*pi)^3*kama);
%%%%%%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%%
K_borzman=1.38*10^(-23);                %������������
T_z=290;                                %�����¶�(K)
Bn=0.8*10^6;                             %��������(Hz)
Fn=4.5;                                 %����ϵ��(db)
noise_parameter.Pn=K_borzman*T_z*Bn*10^(Fn/10); %��������(W)
%%%%%%%%%%%%%%%%%%%�Ӳ��ķ������ģ��%%%%%%%%%%%%%%%%%%%%%%
ampfluct_parameter.model=2;              %�Ӳ��ķ�������ֲ���1ΪLog-normal�ֲ���2ΪWeibull�ֲ� 
ampfluct_parameter.miu=0.2;              %Log-normal�ֲ�����
ampfluct_parameter.sigma=1.2;           %Log-normal�ֲ�����
ampfluct_parameter.b=1.2;                %Weibull�ֲ�����
ampfluct_parameter.c=1.8;                %Weibull�ֲ�����
ampfluct_parameter.u=3;                %K�ֲ�����
ampfluct_parameter.v=7;                %K�ֲ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����״����в���
SBR_parameter.H=800e3;                %����߶�(m)
Re=6378e3;                                 %����뾶
u=398600.8e9;                              %������������
alpha1=0/180*pi;%���µ�γ��
eta_i=45/180*pi;%������
Ve=465.1;%�������ٶ�
Deta=Ve/SBR_parameter.V*(1+SBR_parameter.H/Re);
if isRotation==1
   phi_c=fun_CrabAngle(alpha1/pi*180,eta_i/pi*180,SBR_parameter.H); %������ƫ���� Բ���
   rho_c=fun_CrabMagnitude(alpha1/pi*180,eta_i/pi*180,SBR_parameter.H);    %ƫ������ 
else
    phi_c=0;
    rho_c=1;
end

%% �������
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

%% ���÷�������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

natenna_angle=0*pi/180;       %����б���
crab_theta=0;                         %ƫ����
cita0=pi*(90)/180;                    %���䲨����λ��
phi0=pi*(0)/180;                    %���䲨��������
delta_d_cell=0;                         %��Ԫ���

%%%%%%%%%%%%%%%%%%%�Ӳ��ķ������ģ�Ͳ���%%%%%%%%%%%%%%%%%%%%%%
clutter_amp_fluct_flag=0;               %�Ӳ����������־λ
ampfluct_parameter.model=2;       %�Ӳ��ķ�������ֲ���1ΪLog-normal�ֲ���2ΪWeibull�ֲ�
clutter_freq_fluct_flag=0;              %�Ӳ�Ƶ�������־λ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_noise =1;
%% �����Ӳ�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ru=C/2/fr;                                                          %��ģ������
r_num=floor(ru/delta_r);                                            %��ģ�������ڵľ�������%floor(ru/delta_r)
r_start=1;                                                          %��ʼ������floor(r_min/delta_r)
L = ceil((r_max-H)/ru);                                             %����ģ������
if isRu==0
    L=1;
end
reflect_cell_num=floor(pi/radar_parameter.cita_3db);                %�Ⱦ��뻷���ֵķ��䵥Ԫ��
cita=pi*(0:reflect_cell_num-1)/reflect_cell_num+natenna_angle;      %�Ӳ����䵥Ԫ����λ��
Sys_DOF = arry_col_num*arry_row_num*pulse_num;                      %ϵͳ���ɶ�
% Clutter=zeros(arry_row_num*arry_col_num*pulse_num,r_num+1);       %�Ӳ��洢����(��Ԫ����
% C_Pt=zeros(1,r_num);                                                                            %�Ӳ���ֵ���ʴ洢
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for kk=1:1
Clutter=zeros(arry_row_num*arry_col_num*pulse_num,r_num-r_start+1);       %�Ӳ��洢����(��Ԫ����
C_Pt=zeros(1,r_num-r_start+1);                                                                            %�Ӳ���ֵ���ʴ洢  
R_temp = zeros(arry_row_num*arry_col_num*pulse_num,arry_row_num*arry_col_num*pulse_num);              %�Ӳ�Э���������ʱ�洢����
% for n=0:r_num-r_start%n=1:r_num
for n=0:0%n=1:r_num
    rp=(r_start+n)*delta_r;                                                                                                                                        %�Ӳ����䵥Ԫ�������
    
    Clutter_temp=zeros(arry_row_num*arry_col_num,pulse_num);                                                                        %�Ӳ���ʱ�洢����(��Ԫ����
    Clutter_temp_ch=zeros(arry_col_num,pulse_num);                                                                                           %�Ӳ���ʱ�洢����(ͨ������
    
    
    for l = 0:L-1%L-1
        l
        r=H+rp+ru*l;                                                                                                                                                     %����ģ�����������
        phi=asin(H/r+(r.^2-H^2)/(2*re*r));                                                                                                             %������
        
%         F_t=function_Ft(cita,phi,cita0,phi0);                                                                                                                %���䷽��ͼ
        [F_t(l+1,:),F_r]=Fun_trans_pattern(cita0,phi0,natenna_angle,crab_theta,cita,phi,lambda,natenna_parameter);                                  %����ͽ��շ���ͼ
        P_r=P_r0/r^4;                                                                                                                                              %���书��
        Sigma=Fun_clutter_reflect_model(r,phi,SBR_parameter,radar_parameter,signal_parameter,entironment_parameter);   %�״��Ч�����
        clutter_amp=sqrt(P_r*Sigma).*F_t(l+1,:);                                                                                                                                                %�Ⱦ��뻷�Ӳ��źŵķ���
        
        clutter_amp_fluct=Fun_clutter_amp_fluct(clutter_amp_fluct_flag,reflect_cell_num,ampfluct_parameter,signal_parameter);   %�Ⱦ��뻷����Ԫ���Ӳ��źŵķ������
        clutter_freq_fluct=Fun_clutter_freq_fluct(clutter_freq_fluct_flag,reflect_cell_num,entironment_parameter,SBR_parameter,signal_parameter,radar_parameter);    %�Ӳ��źŵ�Ƶ�����
        
        Phase_rand=exp(1j*pi*(1-2*rand(1,reflect_cell_num)));                                                                                %�Ⱦ��뻷����Ԫ���Ӳ��źŵ������λ
        
        %F_t=ones(1,reflect_cell_num);
        %arry_cell_err=(1+delta_d_cell*randn(arry_row_num,arry_col_num)).*exp(j*pi*delta_d_cell*(unifrnd(-pi,pi,arry_row_num,arry_col_num)));%������Ԫ���
        
        C_0=(clutter_amp.*clutter_amp_fluct.*Phase_rand).'*ones(1,pulse_num).*clutter_freq_fluct;                                     %��׼��Ԫ���յ�pulse_num������Ӳ����ݣ���ʱƵ��
        t_frq=2*rho_c*cos(phi)*cos(cita+phi_c)*V/lambda;         %ʱƵ
        t_ph=exp(1j*2*pi*(t_frq.'*(0:pulse_num-1)/fr));                                                                                            %ʱƵ��λ
        C_z=C_0.*t_ph;                                                                                                                                          %��׼��Ԫ���յ�pulse_num������Ӳ����ݣ���ʱƵ��
        C_Pt(n+1)=clutter_amp*clutter_amp'+C_Pt(n+1);
        
        %% ��ͨ���ϳ�����
        for m=1:reflect_cell_num
            s_frq_azimuth=cos(phi)*cos(cita(m)-natenna_angle)*arry_distance/lambda;                          %��λά��Ƶ
            s_frq_pitch=sin(phi+crab_theta)*arry_distance/lambda;                                %����ά��Ƶ
            s_ph_azimuth=exp(-1j*2*pi*s_frq_azimuth.'*(0:arry_col_num-1));                %��Ƶ��λά��λ
            s_ph_pitch=exp(-1j*2*pi*s_frq_pitch.'*(0:arry_row_num-1));                        %��Ƶ����ά��λ
            arry_cell_err=(1+delta_d_cell*randn(arry_row_num*arry_col_num,1)).*exp(j*pi*delta_d_cell*(unifrnd(-1,1,arry_row_num*arry_col_num,1)));%������Ԫ���
            s_ph(m,:)=kron(s_ph_pitch,s_ph_azimuth).*arry_cell_err.';                                                      %�ϳɵĿ�Ƶ��λ
            %%������ʵ��Э�������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if (r_start+n)==ceil((r_num+r_start)/2)
            X_temp = kron((C_0(m,:).*t_ph(m,:)).',s_ph(m,:).');
            R_temp = X_temp*X_temp' + R_temp;
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        Clutterz=s_ph.'*(C_z);                                                                                  %������Ԫ���յ���ǰ���뻷��pulse_num������Ӳ�����
        Clutter_temp=Clutter_temp+Clutterz;                                                       %�����о���ģ�����Ӳ����ݽ����ۻ�
        
    end
    
    %% ��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flag_noise == 1
        Z =randn(Sys_DOF,1);
        gauss_noise = sqrt(noise_parameter.Pn)*Z.*exp(1j*unifrnd(-pi,pi,Sys_DOF,1));%������˹����
    else
        gauss_noise = zeros(Sys_DOF,1);
    end
    
    %% �������    
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


% %% �Ӳ������׷���
% 
% % ����Э��������ֵD����ͼ
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
figure;plot(D_db','*-');grid on;title('�Ӳ�Э��������������');%-max(D_db)
% % 
% % % % �����Ӳ�������P����ͼ
% % Nfs=180;
% % fs = (-Nfs:1:Nfs)/Nfs;
% % Nfd=180;
% % fd = (-Nfd:1:Nfd)/Nfd;
% % fd=2*cos(phi0)*cos(pi*fs+natenna_angle+crab_theta)*V/lambda/fr;
% % R_inv=inv(R);
% % for p=1:length(fs)
% %     s_ph_azimuth=exp(-j*pi*fs(p)*(0:arry_col_num-1));                         %��Ƶ��λά��λ
% %     s_ph_pitch=exp(-j*pi*fs(p)*(0:arry_row_num-1));                                 %��Ƶ����ά��λ
% %     Ss=kron(s_ph_pitch,s_ph_azimuth);                                                                   %�ϳɵĿ�Ƶ��λ
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
% % title('�Ӳ�������');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB');  %�����Ӳ�������ͼmesh(fd,fs,P_db);cos(\psi)
% % figure;contour(fd,fs,P_db,20);   %�����Ӳ������׵ȸ���ͼ
% % title('�Ӳ������׵ȸ���');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB');  
% %
% %

