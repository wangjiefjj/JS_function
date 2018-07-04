function [ Rk ] = fun_GenerateR( Nt, Nr, Np, Num, CNR, EL, beta, d,gamma, CrabA, CrabM)
%FUN_GENERATER �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% ����Э�������, �����޾���ģ��, ���е�����ת
% Nt: ����������
% Nr�� ����������
%Np�� ������
% Num�� �Ӳ������
% CNR: ����� 
% EL: �����ǣ�rad
% CrabA: ��ƫ��
% CrabM����Ƭ����
if nargin<10    %Ĭ�ϵ�������ת
   CrabA = 0;
   CrabM = 1;
end
D = Nt*Nr*Np;
Rk = zeros(D,D);
AZ_du = linspace(0,179,Num);               %�����Ӳ���
% AZ_du = 90;
AZ = AZ_du/180*pi;
%% ��������
steering_angle = 0; %% ����ָ��
win = chebwin(Num);
for k=1:length(AZ)  
    AF(k) = sum(exp(-1i*pi*d*(0:Nr*Nt-1)*(sin(AZ(k)) ...
                  - sin(steering_angle*pi/180))));
end
AF(1) = AF(2);

AF = diag(abs(AF));

cmj = sin(EL)*cos(AZ);       %����׶��
fspc = 0.5*d*cmj;             %�ռ�Ƶ��
wdc = 0.5*beta*d * CrabM * sin(EL)*cos(AZ+CrabA);   %�Ӳ���1��������
sc = zeros(D,length(cmj));
for i_Az = 1:length(cmj)                   %������          
    a = exp(1j*(0:Nr-1)*2*pi*fspc(i_Az)).';         %���տռ䵼��ʸ�� 
    b = exp(1j*(0:Nt-1)*2*pi*gamma*fspc(i_Az)).';    %����ռ䵼��ʸ��  
    c = exp(1j*(0:Np-1)*2*pi*wdc(i_Az)).';           %ʱ�䵼��ʸ�� 
    sc(:,i_Az) =  kron(c,kron(b,a));             %�Ӳ��ĵ���ʸ��%AF(i_Az) *
end 
Rk = sc*AF*sc';%*AF
Ac=(10^(CNR/10))^0.5; % ����������Ϊ1
Rk=D*Rk/sum(eig(Rk))*Ac^2;
Rk=Rk+eye(D); 
end

