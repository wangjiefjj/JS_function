%% ��С�ɼ���ٶ�
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 9600e6; %450e6   9600e6                 %��Ƶ Hz
C = 299792458;                            %���� m/s
lambda = C/fo;                       %���� m
Nr = 32;                             %������Ԫ����
Nt = 1;                              %������Ԫ����
Np = 32;                             %���������
D = Nt*Nr*Np;                          %���ɶ� 
Re = 6373e3;                          %����뾶m
H = 700e3;  %700e3  %9e3                 %SBR�߶� m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR�ٶ�m/s
end                
d = 1;%                             %��һ����Ԫ���
gamma = 16;                         %������տռ����
PRF = 500;                         %�����ظ�Ƶ��Hz��500~2000
Tr =1/PRF;                         %�����ظ����s
% beta = 19.47;
beta = Vp*2/PRF/(d*lambda/2);      %%���ϵ�� 
alpha1 = 20/180*pi;                %����γ��
eta = 90/180*pi;                   %�������
%% ����
Rs = 130d3;
if H>=500e3
    R = 1301e3;                               %�ؾ�m
    Rs = fun_R2Rs(H, R);
end

%% ����ģ��
% deltaR = fun_deltaR(H, R, Rs, Tr);
% R = [R];  %,R+2*dR,R+3*dR,R+4*dR
%%
EL = asin(H/Rs);         
if H>=500e3
    EL_du = fun_ELAngle(H,R);             %����������
    EL = EL_du/180*pi;
end

%% ƫ���Ƿ���
CrabA = (fun_CrabAngle( alpha1,eta, H)); %ƫ����
CrabM = fun_CrabMagnitude( alpha1,eta, H);%ƫ������
%% �Ӳ�Э����
Num = 360;
CNR = 40;
Rk = zeros(D,D);
Rk_r = zeros(D,D);
for i = 1:length(EL)
    Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,lambda,gamma);%,CrabA,CrabM
    Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma,lambda,CrabA,CrabM);    
end
figure()
mesh(abs(Rk))
iRk = inv(Rk);  
iRk_r = inv(Rk_r); 

Train = fun_TrainData('p',D,2*D,Rk,2,0.5,1);
RSCM = fun_NSCMN(Train);
iRSCM = inv(RSCM);

%% Ŀ��
wd = linspace(-0.5,0.5,200);               %�����1��������
wd_real = wd*PRF/2*lambda/2;                %����ʵ���ٶ�m/s
cmj = wd./beta/(d*lambda^2);
fsp = d*cmj;

for i_wd = 1:length(wd)
    i_wd
    a = exp(1j*(0:Nr-1)*2*pi*fsp(i_wd)).';               %���տռ䵼��ʸ��
    b = exp(1j*(0:Nt-1)*2*pi*gamma*fsp(i_wd)).';         %����ռ䵼��ʸ��
    c = exp(1j*(0:Np-1)*2*pi*wd(i_wd)).';            %ʱ�䵼��ʸ��
    s = kron(c,kron(b,a));                                  %�����ĵ���ʸ��
    SINR_ideal(i_wd) = (abs(s'*iRk*s))^2;          % SINR_ideal ʽ6.37   
    SINR_ideal_r(i_wd) = (abs(s'*iRk_r*s))^2;      % SINR_ideal ʽ6.37   
    SINR(i_wd) = (abs(s'*iRSCM*s))^2/abs(s'*iRSCM*Rk*iRSCM*s);
end
% figure()
% plot(wd,10*log10(SINR_ideal/max(SINR_ideal)))%/max(SINR_ideal)
% hold on
% plot(wd,10*log10(SINR_ideal_r/max(SINR_ideal_r)))%/max(SINR_ideal)
% plot(wd,10*log10(SINR/max(SINR)),'r')%/max(SINR)

