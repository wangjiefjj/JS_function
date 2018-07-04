%% ͼ6.19�� ��λ��90��fsp=0����RDͼ���о���ģ������ת
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9600e6; %450e6   9600e6                 %��Ƶ Hz
C = 299792458;                            %���� m/s
lambda = C/fo;                      %���� m
Nr = 16;                             %������Ԫ����
Nt = 4;                              %������Ԫ����
Np = 16;                             %���������
D = Nt*Nr*Np;                          %���ɶ� 
Re = 6373e3;                          %����뾶m
H = 700e3;  %700e3  %9e3                 %SBR�߶� m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR�ٶ�m/s
end                
d = 60;%                           %��һ����Ԫ���
gamma = 16;                         %������տռ����
PRF = 400;                         %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                         %�����ظ����s
% beta = 19.47;
beta = Vp*2/PRF/(d*lambda/2);      %%���ϵ�� 
alpha1 = 20/180*pi;                  %����γ��
eta = 90/180*pi;                    %�������

%% ƫ���Ƿ���
CrabA = (fun_CrabAngle( alpha1,eta, H)); %ƫ����
CrabM = fun_CrabMagnitude( alpha1,eta, H);%ƫ������

%% Ŀ������
RR = 1000e3:10e3:2400e3;                   %�ؾ�m 
Rs = fun_R2Rs(H,RR);
Graze = fun_GrazeAngle(H, RR, Rs)/180*pi;
dR =  C*Tr/2 *sec(Graze);
%% ģ������͸�����
EL_du(1,:) = fun_ELAngle(H,RR);           %����������
EL(1,:) = EL_du(1,:)/180*pi;
for i = 2:2:4
    RRt = RR + (i-1)*dR;
    EL_du(i,:) = fun_ELAngle(H,RRt);           %����������
    EL(i,:) = EL_du(i,:)/180*pi;
    RRt = RR - (i-1)*dR;
    EL_du(i+1,:) = fun_ELAngle(H,RRt);           %����������
    EL(i+1,:) = EL_du(i,:)/180*pi;   
end
%% �Ӳ�
Num = 180;                          %��λ�ֿ�
CNR = 40;                   %�����dB
wd = linspace(-0.5,0.5,1000);       %�����1��������
%SINR
fsp = 0;
Pt = zeros(length(RR),length(wd));
Ptr = zeros(length(RR),length(wd));
for i = 1:length(RR)
    i
    Rk = zeros(D,D);
    Rk_r = zeros(D,D);
    for j = 1:5
        Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(j,i), beta, d,gamma);%,CrabA,CrabM
        Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(j,i), beta, d,gamma,CrabA,CrabM); 
    end
    for i_wd = 1:length(wd)
        a = exp(1j*(0:Nr-1)*2*pi*fsp).';               %���տռ䵼��ʸ��
        b = exp(1j*(0:Nt-1)*2*pi*gamma*fsp).';         %����ռ䵼��ʸ��
        c = exp(1j*(0:Np-1)*2*pi*wd(i_wd)).';            %ʱ�䵼��ʸ��
        s = kron(c,kron(b,a));                                  %�����ĵ���ʸ��
        Pt(i,i_wd)= (abs(s'*Rk*s));                  % ƥ���˲�SINR 6.27
        Ptr(i,i_wd)= (abs(s'*Rk_r*s));                  % ƥ���˲�SINR 6.27
    end
end

for i_wd = 1:length(wd)
    a = exp(1j*(0:Nr-1)*2*pi*fsp).';               %���տռ䵼��ʸ��
    b = exp(1j*(0:Nt-1)*2*pi*gamma*fsp).';         %����ռ䵼��ʸ��
    c = exp(1j*(0:Np-1)*2*pi*wd(i_wd)).';            %ʱ�䵼��ʸ��
    s = kron(c,kron(b,a));                                  %�����ĵ���ʸ��
    Pt(i,i_wd)= (abs(s'*Rk*s));                  % ƥ���˲�SINR 6.27
    Ptr(i,i_wd)= (abs(s'*Rk_r*s));                  % ƥ���˲�SINR 6.27
end

%% figure
figure(1)
[X,Y]=meshgrid(wd,RR);
mesh(X,Y,10*log10(abs((Pt))))%10*log10
xlabel('��һ��������')
ylabel('����/km')
view(-0,90)
title(['R-Dͼ��H=',num2str(H),'m'])
figure(2)
[X,Y]=meshgrid(wd,RR);
imagesc(wd,RR/1e3,10*log10(abs((Pt))))
xlabel('��һ��������')
ylabel('����/km')
view(-0,-90)
title(['R-Dͼ��������H=',num2str(H),'m'])
figure(3)
imagesc(wd,RR/1e3,10*log10(abs((Ptr))))
xlabel('��һ��������')
ylabel('����/km')
view(-0,-90)
title(['R-Dͼ����ת��������5����H=',num2str(H),'m'])