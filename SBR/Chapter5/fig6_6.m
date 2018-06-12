%% ͼ6.6�� ��λ=89.5��ʹ��2D�����γ�ʱ���Ӳ���,���Դ�����ת�Ͳ�������ת���ܻ�
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1.25e9;                                  %��Ƶ Hz
C = 3e8;                                    %���� m/s
lambda = C/fo;                              %���� m
N = 12;                             %��Ԫ����
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
alpha1 = 45/180*pi;                  %����γ��
eta = 90/180*pi;                    %�������

%% Ŀ������
RR = 400:10:1100;                    %�ؾ�km 
wd = -1:0.001:1;                      %����Ķ����շ�Χ
P = zeros(length(RR),length(wd));
for Az = 89.5:0.1:89.5                   %��λ�� 
    Az
Az = Az/180*pi;
EL = fun_ELAngle(H,RR)/180*pi;       %������
CrabA = (fun_CrabAngle( alpha1,eta, H)); %ƫ����
CrabM = fun_CrabMagnitude( alpha1,eta, H);%ƫ������
% cmj = sin(EL)*cos(Az);                 %����׶��
cmj = CrabM*sin(EL)*cos(Az+CrabA);       %����׶��(������ת)
L = length(cmj);                       %���뵥Ԫ

%% �Ӳ�RDͼ
wdc = beta0*cmj;                 %�Ӳ���1��������
S = zeros(M*N,L);
tic
Pt = zeros(length(cmj),length(wd));
for i = 1:length(cmj)
    b = exp(-1j*(0:M-1)*pi*wdc(i)).';           %ʱ�䵼��ʸ��
    a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %�ռ䵼��ʸ�� 
    S(:,i) = kron(b,a);                         %�Ӳ��ĵ���ʸ��
    X(:,i) = awgn(S(:,i),0);
    R = fun_SCMN(X(:,i));                       %��ͬ�����Ӳ�Э����
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).';         %ʱ�䵼��ʸ��
        a = exp(-1j*(0:N-1)*pi*d*cmj(i)).';         %�ռ䵼��ʸ�� 
        s = kron(b,a);                              %�����ĵ���ʸ��
        Pt(i,i_wd)= (abs(s'*R*s));                  % ƥ���˲�SINR 6.27
    end
end
toc
P = P+Pt;
end
%% figure
figure(1)
[X,Y]=meshgrid(wd,RR);
mesh(X,Y,abs((P)))
xlabel('��һ��������')
ylabel('����/km')
view(-0,90)