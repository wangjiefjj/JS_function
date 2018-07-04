%% �Ƕ�-��������������Ӳ���
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                                  %��Ƶ Hz
C = 3e8;                                    %���� m/s
lambda = C/fo;                              %���� m
N = 16;                             %��Ԫ����
P = 4;
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 7000;                            %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�m/s
d = 63.9;                           %��һ����Ԫ���
% d = 8.3;                           %��һ����Ԫ���
gamma = 5;                              %������տռ����
PRF = 400;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
c = 3e5;                            %����km/m
% beta = 19.47;
beta = Vp*4/PRF/lambda;
beta0 = beta * d;

% beta0 = 2*Vp*Tr/lambda/2;
alpha1 = 45/180*pi;                  %����γ��rad
eta = 90/180*pi;                    %�������rad
Az = linspace(-27/180*pi,27/180*pi,100);
% Az = linspace(0.4,0.42,100);
cosAz = cos(Az);
% cosAz = linspace(-0.05,0.05,100);
%% R=500
RR = 500;                           %�ؾ�km 
EL = fun_ELAngle(H,RR)/180*pi;       %������
wd = -linspace(-1,1,200);                       %��1��������
% Pt = zeros(length(wd),length(Az));
tic
for i_Az = 1:length(cosAz)                 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %����׶��
    %% �Ƕ�-��������������Ӳ��� 
    %% �Ӳ�Э����
    wdc = beta0*cmj;                         %�Ӳ���1��������
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %���տռ䵼��ʸ�� 
    b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).';         %����ռ䵼��ʸ��  
    c = exp(-1j*(0:M-1)*pi*wdc).';           %ʱ�䵼��ʸ�� 
    sc = kron(c,kron(b,a));             %�Ӳ��ĵ���ʸ��
    sc = awgn(sc,20);
    Rk = fun_SCMN(sc);                       %�Ӳ�Э����
    %% ��cos(Az(i_Az))��λ���µ��Ӳ�������
    for i_wd = 1:length(wd)
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %���տռ䵼��ʸ��
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).'./sqrt(P);         %����ռ䵼��ʸ��
        c = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %ʱ�䵼��ʸ��
        s = kron(c,kron(b,a));                           %�����ĵ���ʸ��
        Pt(i_wd,i_Az) = (abs(s'*Rk*s));          % �Ӳ�������6.28    
    end
end
toc
%% figure
figure(1)
imagesc(Az/pi*180,wd,abs((Pt)))
xlabel('\theta_{AZ}/deg')
ylabel('��һ��������')
title(['H=',num2str(H),'km, R=',num2str(RR),'km'])
%% R=1400
RR = 1400;
EL = fun_ELAngle(H,RR)/180*pi;             %������
for i_Az = 1:length(cosAz)                 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);         %����׶��
    %% �Ƕ�-��������������Ӳ��� 
    %% �Ӳ�Э����
    wdc = beta0*cmj;                         %�Ӳ���1��������
    a = exp(-1j*(0:N-1)*pi*d*cmj).';         %���տռ䵼��ʸ�� 
    b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).';         %����ռ䵼��ʸ��  
    c = exp(-1j*(0:M-1)*pi*wdc).';           %ʱ�䵼��ʸ�� 
    sc = kron(c,kron(b,a));             %�Ӳ��ĵ���ʸ��
    sc = awgn(sc,20);
    Rk = fun_SCMN(sc);                       %�Ӳ�Э����
    %% ��cos(Az(i_Az))��λ���µ��Ӳ�������
    for i_wd = 1:length(wd)
        a = exp(-1j*(0:N-1)*pi*d*cmj).'./sqrt(N);         %���տռ䵼��ʸ��
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj).'./sqrt(P);         %����ռ䵼��ʸ��
        c = exp(-1j*(0:M-1)*pi*wd(i_wd)).'./sqrt(M);           %ʱ�䵼��ʸ��
        s = kron(c,kron(b,a));                           %�����ĵ���ʸ��
        Pt2(i_wd,i_Az) = (abs(s'*Rk*s));          % �Ӳ�������6.28    
    end
end

%% figure
figure(2)
imagesc(Az/pi*180,wd,abs((Pt2)))
xlabel('\theta_{AZ}/deg')
ylabel('��һ��������')
title(['H=',num2str(H),'km, R=',num2str(RR),'km'])