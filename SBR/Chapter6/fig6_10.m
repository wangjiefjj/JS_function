%% ����ͼ
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.6e9;                         %��Ƶ Hz
C = 3e8;                            %���� m/s
lambda = C/fo;                      %���� m
N = 16;                             %��Ԫ����
P = 4;
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 1250;                           %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�m/s
d = 1;                            %��һ����Ԫ���
% d = 8.3;                          %��һ����Ԫ���
gamma = 5;                         %������տռ����
PRF = 400;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
c = 3e5;                            %����km/m
% beta = 19.47;
beta = Vp*4/PRF/lambda;
beta0 = beta * d;
%% ����ͼ 
theta0_AZ = 10/180*pi;                 %��λ��ָ��
% theta0_EL = 10/180*pi;                 %������ָ��
theta_AZ_du = (-100:1:100);
theta_AZ=(theta_AZ_du)./180*pi;           %��λ��������Χ�Ͳ���
theta_EL=(0:40)./180*pi;                %������������Χ�Ͳ���
% w = ones(N*P,1);
% window = taylorwin(N*P);
% s = s.*window;
for i_EL = 1:length(theta_EL)
    i_EL
    c_s = cos(theta_EL(i_EL))*sin(theta0_AZ);   %�ź�׶�� sin(theta0_EL)*
    a = exp(-1j*(0:N-1)*pi*d*c_s).';   %
    b = exp(-1j*(0:P-1)*pi*gamma*d*c_s).';
    w = kron(b,a);                  %Ȩֵ
    w = w/sqrt(w'*w);               %��һ��
    cmj = cos(theta_EL(i_EL))*sin(theta_AZ); 
    for i_AZ = 1:length(theta_AZ)
        a = exp(-1j*(0:N-1)*pi*d*cmj(i_AZ)).';               %���տռ䵼��ʸ�� 
        b = exp(-1j*(0:P-1)*pi*gamma*d*cmj(i_AZ)).';         %����ռ䵼��ʸ��
        s = kron(b,a);                                       %�ź�
        s = s/sqrt(s'*s);                                 %��һ��
%         w = w.*window;
        A(i_AZ,i_EL) = w'*s;                          %����ͼ(ÿһ����һ����λ��)       
    end
end
figure()
plot(theta_AZ_du,10*log10(abs(A(:,1))))
[hang,lie] = size(A);
if lie>1
    figure()
    mesh(10*log10(abs(A)))
end
% axis([-100,100,-1e-7,0])