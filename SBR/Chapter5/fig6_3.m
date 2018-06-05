%% ͼ6.3������500km�����������Ķ�����Ƶ�ʣ��е�����ת
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;                             %��Ԫ����
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 506;                            %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�km/s
d = 13.4;                           %��һ����Ԫ���
PRF = 2000;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
c = 3e5;                            %����km/m
beta = 19.47;
beta0 = beta * d;
alpha1 = 30/180*pi;                  %����γ��
eta = 45/180*pi;                    %�������
%% Ŀ��
R=500;                              %��ʵĿ�����km
Rs = fun_R2Rs(H,R);
Graze = fun_GrazeAngle(H, R, Rs)/180*pi;
dR =  c*Tr/2 *sec(Graze);   
RR=R+(-2:4)*dR;                     %ģ������ 7��
%% ���������
c = -1:0.01:1;                       %����׶��
Az = (0:180)/180*pi;
ELm = fun_ELAngle(H,RR)/180*pi;
% wd=0;
for i = 1:length(ELm)
    cosAz = c./sin(ELm(i));
    wd(i,:) = beta0*sin(ELm(i)).*cosAz;
end
CrabA = (fun_CrabAngle( alpha1,eta, H));
CrabM = fun_CrabMagnitude( alpha1,eta, H);
% wdr=0;
for i = 1:length(ELm)
    Az = (acos(c./sin(ELm(i))));
    wdr(i,:) = beta0*CrabM*sin(ELm(i)).*cos(Az+CrabA);
end
wdr = real(wdr)+imag(wdr);
%% figure
figure(1)
hold on
plot(c,wd','b')
plot(c,(wdr.'),'r')
