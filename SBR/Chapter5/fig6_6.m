%% ͼ6.6���޵�����ת�� ��λ=89.5��ʹ��2D�����γ�ʱ���Ӳ���
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 32;                             %��Ԫ����
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 506;                            %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�
d = 13.4;                           %��һ����Ԫ���
PRF = 500;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
c = 3e5;                            %����km/m
beta = 19.47;
beta0 = beta * d;
% alpha1 = 30/180*pi;                  %����γ��
% eta = 45/180*pi;                    %�������
%% Ŀ������
RR = 400:10:1100;                    %�ؾ�km 
Az = 89.5/180*pi;                    %��λ�� 
EL = fun_ELAngle(H,RR)/180*pi;       %������
c = sin(EL)*cos(Az);                 %����׶��
L = length(c);                       %���뵥Ԫ

% c = -1:0.01:1;
% Az = 89.5/180*pi;                      %��λ�� 
% L = length(c);                       %���뵥Ԫ
% RR = linspace(400,1100,L);   

%% Ŀ�굼��ʸ������
wdc = beta0*c/PRF;                 %��1��������
S = zeros(M*N,L);
% w = hamming(L);
%�Ӳ�����
v = 1;
for i = 1:length(c)
    b = exp(-1j*(0:M-1)*pi*wdc(i)).';         %ʱ�䵼��ʸ��
    a = exp(-1j*(0:N-1)*pi*d*c(i)).';         %�ռ䵼��ʸ�� 
    S(:,i) = gamrnd(v,1/v,1,1)*kron(b,a);     %�Ӳ��ĵ���ʸ��
%     X(:,i) = awgn(S(:,i),0);
end
R = fun_SCMN(S);
invR = inv(R);
wd = -1:0.01:1;
for i_c = 1:length(c)
    for i_wd = 1:length(wd)
        b = exp(-1j*(0:M-1)*pi*wd(i_wd)).';         %ʱ�䵼��ʸ��
        a = exp(-1j*(0:N-1)*pi*d*c(i_c)).';        %�ռ䵼��ʸ�� 
        s = kron(b,a);
        P(i_c,i_wd)= s'*R*s;
    end
end
%% figure
figure(1)
[X,Y]=meshgrid(wd,RR);
mesh(X,Y,abs(fftshift(P)))
xlabel('��һ��������')
ylabel('����/km')
