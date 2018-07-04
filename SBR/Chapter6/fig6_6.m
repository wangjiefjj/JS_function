%% ͼ6.6�� ��λ=89.5��ʹ��2D�����γ�ʱ���Ӳ���,���Դ�����ת�Ͳ�������ת���ܻ�
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                                  %��Ƶ Hz
C = 3e8;                                    %���� m/s
lambda = C/fo;                              %���� m
N = 16;                             %��Ԫ����
P = 4;
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 700e3;                            %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�m/s
d = 1;                           %��һ����Ԫ���
% d = 8.3;                           %��һ����Ԫ���
gamma = 5;                              %������տռ����
PRF = 400;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
c = 3e5;                            %����km/m
% beta = 19.47;
beta = Vp*2/PRF/(d*2*lambda);
beta0 = beta * d;

alpha1 = 45/180*pi;                  %����γ��
eta = 90/180*pi;                    %�������

%% Ŀ������
RR = 800e3:10e3:2400e3;                    %�ؾ�km 
wd = -1:0.01:1;                      %����Ķ����շ�Χ
PP = zeros(length(RR),length(wd));
for Az = 0:0.1:0                   %��λ��/deg
    Az
Az = Az/180*pi;
EL = fun_ELAngle(H,RR)/180*pi;       %������
CrabA = (fun_CrabAngle( alpha1,eta, H)); %ƫ����
CrabM = fun_CrabMagnitude( alpha1,eta, H);%ƫ������
cmj = sin(EL)*cos(Az);                     %����׶��
% cmj = CrabM*sin(EL)*cos(Az+CrabA);       %����׶��(������ת)
L = length(cmj);                       %���뵥Ԫ

%% �Ӳ�RDͼ
wdc = beta0*cmj;                 %�Ӳ���1��������
S = zeros(P*M*N,L);
tic
Pt = zeros(length(cmj),length(wd));
for i = 1:length(cmj)
    i
    a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %���տռ䵼��ʸ�� 
    b = exp(-1j*(0:P-1)*pi*gamma*d*cmj(i)).';         %����ռ䵼��ʸ��  
    c = exp(-1j*(0:M-1)*pi*wdc(i)).';           %ʱ�䵼��ʸ�� 
    S(:,i) = kron(c,kron(b,a));                         %�Ӳ��ĵ���ʸ��
    X(:,i) = awgn(S(:,i),20);
    R = fun_SCMN(X(:,i));                       %��ͬ�����Ӳ�Э����
    for i_wd = 1:length(wd)
        a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %���տռ䵼��ʸ�� 
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj(i)).';         %����ռ䵼��ʸ��  
        c = exp(-1j*(0:M-1)*pi*wd(i_wd)).';         %ʱ�䵼��ʸ��    
        s = kron(c,kron(b,a));                              %�����ĵ���ʸ��
        Pt(i,i_wd)= (abs(s'*R*s));                  % ƥ���˲�SINR 6.27
    end
end
toc
PP = PP+Pt;
end
%% figure
% figure(1)
% [X,Y]=meshgrid(wd,RR);
% mesh(X,Y,abs((P)))
% xlabel('��һ��������')
% ylabel('����/km')
% view(-0,90)
% title(['R-Dͼ��H=',num2str(H),'km'])
figure(2)
[X,Y]=meshgrid(wd,RR);
imagesc(wd,RR,abs((PP)))
xlabel('��һ��������')
ylabel('����/km')
view(-0,-90)
title(['R-Dͼ��H=',num2str(H),'km'])