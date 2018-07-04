%% ���ղ����γ�
%% ����ͼ
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                         %��Ƶ Hz
C = 3e8;                            %���� m/s
lambda = C/fo;                      %���� m
Nrx = 16;                           %���շ�λ����Ԫ����
Nry = 2;                            %���ո�������Ԫ����
P = 4;                              %����ͨ��
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 1250;                           %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�m/s
drx = 63.9;                         %����x���һ����Ԫ���
dry = 16;                           %����y���һ����Ԫ���
% d = 8.3;                          %��һ����Ԫ���
gamma = 16;                          %������տռ����(x����)
PRF = 400;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
c = 3e5;                            %����km/m
% beta = 19.47;
% beta = Vp*4/PRF/lambda;
% beta0 = beta * d;
%% �����ź�
theta0_AZ = 50/180*pi;                 %��λ��ָ��
theta0_EL = 20/180*pi;                 %������ָ��

theta_AZ_du = (-100:1:100);
theta_EL_du = (-100:1:100);
theta_AZ=(theta_AZ_du)./180*pi;         %��λ��������Χ�Ͳ���
theta_EL=(theta_EL_du)./180*pi;                %������������Χ�Ͳ���
%% Ȩֵ
c_sx = sin(theta0_AZ);   %�ź�׶�� sin(theta0_EL)*
c_sy = sin(theta0_AZ)*sin(theta0_EL);   %�ź�׶�� sin(theta0_EL)*
ax = exp(-1j*(0:Nrx-1)*pi*drx*c_sx).';   %
ay = exp(-1j*(0:Nry-1)*pi*dry*c_sy).';   %
a = kron(ax,ay);
w = a;
% b = exp(-1j*(0:P-1)*pi*gamma*drx*c_sx).';
% w = kron(b,a);                  %Ȩֵ
% w = w/sqrt(w'*w);               %��һ��
%% ����ͼ
for i_EL = 1:length(theta_EL)
    i_EL
    cx =  sin(theta_AZ);
    cy =  sin(theta_AZ)*sin(theta_EL(i_EL)); %sin(theta_EL(i_EL))* sin(theta_AZ)
    for i_AZ = 1:length(theta_AZ)        
        ax = exp(-1j*(0:Nrx-1)*pi*drx*cx(i_AZ)).';               %���տռ䵼��ʸ��x�� 
        ay = exp(-1j*(0:Nry-1)*pi*dry*cy(i_AZ)).';   %
        a = kron(ax,ay);
        s = a;
%         b = exp(-1j*(0:P-1)*pi*gamma*drx*cx(i_AZ)).';         %����ռ䵼��ʸ��
%         s = kron(b,a);                                       %�ź�
%         s = s/sqrt(s'*s);                                 %��һ��
%         w = w.*window;
        A(i_AZ,i_EL) = w'*s;                          %����ͼ(ÿһ����һ����λ��)       
    end
end
figure()
plot(theta_AZ_du,10*log10(abs(A(:,1))))
figure()
plot(theta_EL_du,10*log10(abs(A(1,:))))
[hang,lie] = size(A);
if lie>1
    figure()
    mesh(10*log10(abs(A)))
end
% axis([-100,100,-1e-7,0])