%% R=500km�ĽǶ�-��������������Ӳ���
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9.26e9;                                  %��Ƶ Hz
C = 3e8;                                    %���� m/s
lambda = C/fo;                              %���� m
N = 32;                             %��Ԫ����
M = 16;                             %���������
Re = 6373;                          %����뾶km
H = 1250;                            %SBR�߶�km
Vp = fun_Vp(H);                     %SBR�ٶ�m/s
d = 63.9;                           %��һ����Ԫ���
% d = 10;                            %��һ����Ԫ���
PRF = 500;                          %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                          %�����ظ����s
beta =  Vp*2/PRF/lambda;             %
beta0 = beta * d;

alpha1 = 45/180*pi;                  %����γ��rad
eta = 90/180*pi;                     %�������rad

%% �Ӳ�Э����
RR = 500;                           %�ؾ�km 
EL = fun_ELAngle(H,RR)/180*pi;       %������
Rk = fun_GenerateR(1, N, M, 256, 40, EL, beta0, d);
% Az = linspace(0,180,256)/180*pi;     %��λ��deg
% % cosAz = cos(Az);
% %% ׶�Ǻ��Ӳ�������
% CrabA = (fun_CrabAngle( alpha1,eta, H)); %ƫ����
% CrabM = fun_CrabMagnitude( alpha1,eta, H);%ƫ������
% % cmj = sin(EL).*cos(AZ);         %����׶��
% cosAz = cos(Az+CrabA);



% tic
% %% �Ӳ�Э����
% Rk = zeros(M*N,M*N);
% CNR = 40;                                  %�����
% for i_Az = 1:length(cosAz)                 %������
%     i_Az
%     cmj = sin(EL)*CrabM*cosAz(i_Az);               %����׶��
%     wdc = beta0*cmj;                         %�Ӳ���1��������
%     b = exp(-1j*(0:M-1)*pi*wdc).';           %ʱ�䵼��ʸ��
%     a = exp(-1j*(0:N-1)*pi*d*cmj).';         %�ռ䵼��ʸ�� 
%     sc = kron(b,a);                          %�Ӳ��ĵ���ʸ��
%     Rk = Rk + fun_SCMN(sc);                  %�Ӳ�Э����
% end
% Rk = Rk/max(max(abs(Rk)));
% Rk = Rk*10^(CNR/10) + eye(size(Rk));
% iRk = inv(Rk);

%% �Ƕ�-��������������Ӳ���
% Az = linspace(0,180,256)/180*pi;     %��λ��deg
% cosAz = cos(Az);
cosAz = linspace(-1,1,256);
wd = -linspace(-1,1,200);                       %��1��������
for i_Az = 1:length(cosAz) 
    i_Az
    cmj = sin(EL)*cosAz(i_Az);               %����׶��
    for i_wd = 1:length(wd)
         b = exp(-1j*(0:M-1)*pi*wd(i_wd)).';           %ʱ�䵼��ʸ��
         a = exp(-1j*(0:N-1)*pi*d*cmj).';         %�ռ䵼��ʸ��
         s = kron(b,a);                           %�����ĵ���ʸ��
         Pt(i_wd,i_Az) =(abs(s'*Rk*s));          % �Ӳ�������6.28    
    end
end
%% figure
figure(1)
[X,Y]=meshgrid(cosAz,wd);
mesh(X,Y,((Pt)))
xlabel('cos(\theta)')
ylabel('��һ��������')
view(-0,90)

figure(2)
imagesc(cosAz,wd,((Pt)))
xlabel('cos(\theta)')
ylabel('��һ��������')