%% ȫ���ɶ�stap
clc;clear;close all
isRu=0;
%% ���ƽ̨��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo = 1250e6; %450e6   9600e6                 %��Ƶ Hz
c = 299792458;                            %���� m/s
lambda = c/fo;                      %���� m
Nr = 32;                             %������Ԫ����
Nt = 1;                              %������Ԫ����
Np = 16;                             %���������
D = Nt*Nr*Np;                          %���ɶ� 
Re = 6373e3;                          %����뾶m
H = 850e3;  %700e3  %9e3                 %SBR�߶� m
Vp = 50;
if H >=500e3
    Vp = fun_Vp(H);                      %SBR�ٶ�m/s
end                
d = lambda/2;%                     %��Ԫ���
gamma = 8;                        %������տռ����
fr = 500;                         %�����ظ�Ƶ��Hz��500~2000
Tr =1/fr;                         %�����ظ����s
% beta = 19.47;
beta = Vp*2/fr/d;      %%���ϵ�� 
alpha1 = 20/180*pi;                  %����γ��
eta = 90/180*pi;                    %�������

%% ƫ���Ƿ���
CrabA = (fun_CrabAngle( alpha1,eta, H)); %ƫ����
CrabM = fun_CrabMagnitude( alpha1,eta, H);%ƫ������

%% Ŀ������
R = 1000e3;                                %�ؾ�m 
Rs = fun_R2Rs(H,R);
El0 = fun_ELAngle(H,R);
% %% ģ������͸�����
Rsmax = fun_Rsmax(H);                       %���̽��б��m
Rsu=c/(2*fr);                               %�����ģ��б��m
%% �Ӳ�Э����
Num = 180;                          %��λ�ֿ�
CNR = 40;                          %�����dB
wd = linspace(-0.5,0.5,100);       %�����1��������
cmj = wd./beta/(d*2*lambda);
fspc = d*cmj;
Pt = zeros(length(R),length(wd));
Ptr = zeros(length(R),length(wd));
%����ģ��
Nk1 = floor((Rsmax-Rs)/Rsu);                %ǰ��ģ����
Nk2 = floor((Rs-H)/Rsu);                    %����ģ����
Rsk1 = Rs + (0:Nk1)*Rsu;                    %ǰ��ģ��б��m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %����ģ��б��m
if isRu ==1
    Rsk = [Rsk2,Rsk1];%Rs(i);%                 %ģ��б��m
else
    Rsk = Rs;
end
RRk = fun_Rs2R(H,Rsk);                       %ģ���ؾ�m
elk = fun_ELAngle(H,RRk)/180*pi;             %��ģ�����븩����rad
Rk = zeros(D,D);
Rk_r = zeros(D,D);
for j = 1:length(elk)
    Rk = Rk + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(j), beta, d,gamma, lambda);%,CrabA,CrabM
    Rk_r = Rk_r + fun_GenerateR(Nt, Nr, Np, Num, CNR, elk(j), beta, d,gamma,lambda,CrabA,CrabM); 
end
%% Target Space-Time Steering Vector
azt = 0; elt = El0;                                      % Target azimuth and elevation angles in degrees.
fdt = 100;                                               % Target Doppler Frequency.
omegact = fdt/fr;                                        % Normalized Target Frequency.
fspt = d/lambda*cos(elt*pi/180)*cos(azt*pi/180);     % Target Spatial Frequency.
at = exp(1i*2*pi*fspt*(0:Nr-1));                          % Target Spatial Steering Vector.
bt = exp(1i*2*pi*omegact*(0:Np-1));                       % Target Doppler Steering Vector
vt = kron(bt,at).';                                      % Target Space-Time Steering Vector.

%% Optimum, Fully Adaptive STAP Solution
w = Rk\vt;                                               % Eq. (104)
w_r = Rk_r\vt;
%% Adapted Patterns
az = -90:.5:90;     Laz = length(az);
fd = -150:.5:150;  Lfd = length(fd);
fsp = d/lambda*cos(El0)*cos(az*pi/180);
omega = fd/fr;
Pw1 = zeros(Lfd,Laz);
for m=1:Laz
    for n=1:Lfd
        a = exp(1i*2*pi*fsp(m)*(0:Nr-1));                % Dummy Spatial Steering Vector.(Dummy����)
        b = exp(1i*2*pi*omega(n)*(0:Np-1));              % Dummy Doppler Steering Vector
        v = kron(b,a).';
        Pw1(n,m) = abs(w'*v)^2;
    end
end

%% Normalization:
max_value = max(max(Pw1));
Pw = Pw1/max_value;
%% Plot the Adapted Pattern
figure('NumberTitle', 'off','Name', ...
       'Figure 23a. Example Scenario: Adapted Pattern for Optimum Fully Adaptive STAP', ...
       'Position',[1 1 700 600]);
[Az Doppler] = meshgrid(sin(az*pi/180),fd);
colormap jet;
mesh(Az, Doppler, 10*log10(abs(Pw)));
shading interp;
xlim([-1 1])
ylim([-150 150]);
xlabel('sin(Azimuth)');
ylabel('Doppler Frequency (Hz)');
h = colorbar;
% h = colorbar('YTickLabel',{-80:10:0});
set(get(h,'YLabel'),'String','Relative Power (dB)');