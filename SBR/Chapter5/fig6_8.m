%% R=500km�ĽǶ�-��������������Ӳ���
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1.25e9;                                  %��Ƶ Hz
C = 3e8;                                    %���� m/s
lambda = C/fo;                              %���� m
N = 32;                             %��Ԫ����
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 506;                            %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�m/s
d = 13.4;                           %��һ����Ԫ���
PRF = 500;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
c = 3e5;                            %����km/m
beta = 19.47;
beta0 = beta * d;
% beta0 = 2*Vp*Tr/lambda/2;
alpha1 = 45/180*pi;                  %����γ��rad
eta = 90/180*pi;                    %�������rad
cosAz = linspace(-0.05,0.05,100);
RR = 500;                           %�ؾ�km 
EL = fun_ELAngle(H,RR)/180*pi;       %������
wd = -linspace(-1,1,200);                       %��1��������
% Pt = zeros(length(wd),length(Az));
tic
for i_Az = 1:length(cosAz)                 %������
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %����׶��
    %% �Ƕ�-��������������Ӳ��� 
    %% �Ӳ�Э����
    wdc = beta0*cmj;                         %�Ӳ���1��������
    b = exp(-1j*(0:M-1)*pi*wdc).';           %ʱ�䵼��ʸ��
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %�ռ䵼��ʸ�� 
    sc = kron(b,a);                          %�Ӳ��ĵ���ʸ��
    Rk = fun_SCMN(sc);                       %�Ӳ�Э����
%     Rk = Rk + eye(M*N);
%     iRk = inv(Rk);
    %% ��cos(Az(i_Az))��λ���µ��Ӳ�������
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %ʱ�䵼��ʸ��
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %�ռ䵼��ʸ��
        s = kron(b,a);                           %�����ĵ���ʸ��
        Pt(i_wd,i_Az) = (abs(s'*Rk*s));          % �Ӳ�������6.28    
    end
end
toc
%% figure
figure(1)
imagesc(cosAz,wd,abs((Pt)))
xlabel('cos(\theta)')
ylabel('��һ��������')
title('R=500km')
%% R=1400
RR = 1400;
EL = fun_ELAngle(H,RR)/180*pi;             %������
for i_Az = 1:length(cosAz)                 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %����׶��
    %% �Ƕ�-��������������Ӳ��� 
    %% �Ӳ�Э����
    wdc = beta0*cmj;                         %�Ӳ���1��������
    b = exp(-1j*(0:M-1)*pi*wdc).';           %ʱ�䵼��ʸ��
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %�ռ䵼��ʸ�� 
    sc = kron(b,a);                          %�Ӳ��ĵ���ʸ��
    sc = awgn(sc,40);
    Rk = fun_SCMN(sc);                       %�Ӳ�Э����
    %% ��cos(Az(i_Az))��λ���µ��Ӳ�������
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %ʱ�䵼��ʸ��
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %�ռ䵼��ʸ��
        s = kron(b,a);                           %�����ĵ���ʸ��
        Pt2(i_wd,i_Az) = (abs(s'*Rk*s));          % �Ӳ�������6.28    
    end
end

%% figure
figure(2)
imagesc(cosAz,wd,abs((Pt2)))
xlabel('cos(\theta)')
ylabel('��һ��������')
title('R=1400km')