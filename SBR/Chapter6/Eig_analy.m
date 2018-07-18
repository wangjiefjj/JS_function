%% ����ֵ
clc;clear;close all
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1250e6; %450e6   9600e6 1250e6                 %��Ƶ Hz
C = 3e8;%299792458;                            %���� m/s
lambda = C/fo;                       %���� m
Nr = 16;                             %������Ԫ����
Nt = 1;                              %������Ԫ����
Np = 8;                             %���������
D = Nt*Nr*Np;                          %���ɶ� 
Re = 6373e3;                          %����뾶m
H = 800e3;  %700e3  %9e3                 %SBR�߶� m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR�ٶ�m/s
end                
d = lambda/2;%                             %��һ����Ԫ���
gamma = 16;                         %������տռ����
fr = 500;                         %�����ظ�Ƶ��Hz��500~2000
Tr =1/fr;                         %�����ظ����s
% beta = 19.47;
beta = Vp*2/fr/d;           %%���ϵ�� 
alpha1 = 30;                %����γ��
eta = 70;                   %�������
%% ����
Rs = 130e3;
if H>=500e3
    R =  1.149762283480509e+06;                              %�ؾ�m
    Rs = fun_R2Rs(H, R);
end

%% ����ģ��
Rsmax = fun_Rsmax(H);                       %���̽��б��m
Rsu=C/(2*fr);                               %�����ģ��б��m
Nk1 = floor((Rsmax-Rs)/Rsu);                %ǰ��ģ����
Nk2 = floor((Rs-H)/Rsu);                    %����ģ����
Rsk1 = Rs + (0:Nk1)*Rsu;                    %ǰ��ģ��б��m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %����ģ��б��m
Rsk = [Rsk2,Rsk1];%Rs;%                     %ģ��б��m
Rk = fun_Rs2R(H,Rsk);                       %ģ���ؾ�m
elk = fun_ELAngle(H,Rk)/180*pi;                    %��ģ�����븩����rad
grazek =  fun_GrazeAngle(H,Rk,Rsk)/180*pi;         %��ģ�����������rad
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
Num = 180;
CNR = 40;
Rk = zeros(D,D);
Rk_r = zeros(D,D);
%�޾���ģ��
for i = 1:length(EL)%%EL
    i
    Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma,lambda);%,CrabA,CrabM 
    Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, EL(i), beta, d,gamma,lambda,CrabA,CrabM);
end
E=abs(eig(Rk));
E = 10*log10(sort(E,'descend')).';
figure
hold on
plot(E,'r')
E_r=abs(eig(Rk_r));
E_r = 10*log10(sort(E_r,'descend')).';
plot(E_r,'g')
%�о���ģ��
Rku = zeros(D,D);
Rku_r = zeros(D,D);
for i = 1:length(elk)%%EL
    i
    Rku_r = Rku_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(i), beta, d,gamma,lambda,CrabA,CrabM); 
    Rku = Rku + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(i), beta, d,gamma,lambda);%,CrabA,CrabM    
end
Eu=abs(eig(Rku));
Eu = 10*log10(sort(Eu,'descend')).';
hold on
plot(Eu,'b')
Eu_r=abs(eig(Rku_r));
Eu_r = 10*log10(sort(Eu_r,'descend')).';
plot(Eu_r,'k')
grid on 
legend('�޾���ģ������ת','�޾���ģ������ת','�о���ģ������ת','�о���ģ������ת')