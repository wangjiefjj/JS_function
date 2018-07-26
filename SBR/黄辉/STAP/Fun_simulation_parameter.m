%function Fun_simulation_parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%��������%%%%%%%%%%%%%%%%%%%%%%
entironment_parameter.C=3*10^8;               %���٣�m/s)
entironment_parameter.re=8490000;                %�������ʰ뾶��m)
entironment_parameter.Vm=6;                   %���٣�m/s)
entironment_parameter.scene_parameter=5;      %����������1 ���Ӳ� 2 ɳĮ 3 ũ�� 4 ���� 5 ��ɽ 
entironment_parameter.ss=2;                   %���������1 2 3 4 5������
%%%%%%%%%%%%%%%%%%%%%%�ػ�����%%%%%%%%%%%%%%%%%%%%%%%
airborne_parameter.H=6000;                %�ػ��߶�(m)
airborne_parameter.V=350;                 %�ػ��ٶȣ�m/s)
%%%%%%%%%%%%%%%%%%%%%%�źŲ���%%%%%%%%%%%%%%%%%%%%%%%
signal_parameter.fz=1.5;           %��Ƶ��GHz��
signal_parameter.lambda=0.3/signal_parameter.fz;              %����(m)
signal_parameter.Bs=0.8*10^6;              %�źŴ���
signal_parameter.minu=40*10^(-6);         %����s��
signal_parameter.fs=1*10^6;              %������(��/s)
signal_parameter.fr=4*airborne_parameter.V/signal_parameter.lambda;            %�����ظ�Ƶ�ʣ�Hz),8*airborne_parameter.V/signal_parameter.lambda
signal_parameter.pulse_num=16;           %������
signal_parameter.D_zk=signal_parameter.minu*signal_parameter.fr;             %ռ�ձ�
signal_parameter.D=signal_parameter.Bs*signal_parameter.minu;             %��ѹ��
signal_parameter.delta_r=entironment_parameter.C/signal_parameter.fs/2;        %����ֱ浥Ԫ��m��
%signal_parameter.delta_r=150;           %����ֱ浥Ԫ(km)
%%%%%%%%%%%%%%%%%%%%%���в���%%%%%%%%%%%%%%%%%%%%%%%%
natenna_parameter.arry_distance=0.5*signal_parameter.lambda;  %��Ԫ���
%natenna_parameter.arry_distance=0.5;  %��Ԫ���
natenna_parameter.arry_col_num=16;                     %��Ԫ����
natenna_parameter.Wq_azimuth_db=35;          %�м�Ȩ��db��
natenna_parameter.arry_row_num=1;                     %��Ԫ����
natenna_parameter.Wq_pitch_db=35;           %�м�Ȩ��db)

%%%%%%%%%%%%%%%%%%%%%%�״����%%%%%%%%%%%%%%%%%%%%%%%
Pt=4.8*10^6;                            %�״﷢��ķ�ֵ����(W)
Gt0=1;                                  %�״﷢����������
Gr0=1;                                  %�״������������
kama=10^(8/10);                                 %ϵͳ�������
radar_parameter.r_min=6e3;               %�״��ʼ���þ���(Km)
radar_parameter.r_max=sqrt(2*entironment_parameter.re*airborne_parameter.H);     %�״�ֱ�Ӿ��루m)
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
ampfluct_parameter.model=3;              %�Ӳ��ķ�������ֲ���1ΪLog-normal�ֲ���2ΪWeibull�ֲ� , 3ΪK�ֲ�,K �ֲ�ʱ��U��V(Ҫ��������) 
ampfluct_parameter.miu=0.2;              %Log-normal�ֲ�����
ampfluct_parameter.sigma=1.2;           %Log-normal�ֲ�����
ampfluct_parameter.b=1.2;                %Weibull�ֲ�����
ampfluct_parameter.c=1.8;                %Weibull�ֲ�����
ampfluct_parameter.u=3;                %K�ֲ�����
ampfluct_parameter.v=7;                %K�ֲ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








