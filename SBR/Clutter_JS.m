%% ���� ����״��ʱ����Ӧ�Ӳ����Ƽ���,������ �����µ�
clc
clear
close all
%% ����
Re = 6373e3;                    %Բ����뾶m
c  = 3e8;                 %����m/s 299792458
%% �״�ϵͳ����
fo = 1.25e9;                    %��ƵHz
lambda = c/fo;                  %����
Nr_az = 16;                     %���շ�λ�������
Nr_el = 25;                     %���շ�λ�������
Nt_az = 16;                     %���䷽λ�������
Nt_el = 25;                     %���䷽λ�������
Ti = 14.2e-3;                   %פ��ʱ��s
fr = 5000;                      %�����ظ�����Hz
Np = fr*Ti;                     %���������� 
Pt = 30e3;                      %��ֵ����W
Tp = 20e-6;                     %����s
Pa = Pt*Tp/(1/fr);              %ƽ�����书��
Ls = 6;                         %ϵͳ���dB
B = 1e6;                        %����Hz
D = Tp*B;                       %��ѹ��
Gt = 23;                        %��������dB
Gr = 23;                        %��������dB
d = lambda/2;                   %��Ԫ���m(������)
d_saz = 0.12*25;                %��λ��������m
%% �״�
alpha_r = 30;                   %�״�γ��N,deg
beta_r = 110;                   %�״ﾭ��E,deg
H = 850e3;                      %�߶�m
Vp = fun_Vp(H);                 %ƽ̨�ٶ�
phi = 70;                       %������ 
alpha1 = 110;                   %�״ﾭ��E,deg
beta1 =  30;                    %�״�γ��N,deg
%% Ŀ��
graze = 30;                     %�����deg
alpha2 = 120;                  %Ŀ�꾭��E,deg
beta2  = 25;                   %Ŀ��γ��N,deg
RCS = 5;                        %Ŀ��RCS��dBsm����dB*m^2��
R = ( pi/2 - graze/180*pi - asin( Re/(Re+H) * cos(graze/180*pi) ) ) * Re; %(2.26)����
Rs = fun_R2Rs(H,R);
% ��λ��
Nc = 359;                           %�Ӳ������ 
Az = linspace(0,179,Nc);            %��λ�Ƕ�
El = linspace(0,179,Nc);            %������Ƕ�
% Az = Az/180*pi;
% El = El/180*pi;
EL =  fun_ELAngle(H,R);              %��б�ǣ����߷��߸����ǣ�
LAz = length(Az);
AFt = ones(1,LAz);                  % ��������������
AFr = ones(1,LAz);                  % ��������������

%% �ն�SNR����
k = 1.3806488e-23;                  % ������������ J/K.
T0 = 300;                           % ��׼���������� Kelvin.
Fn = 3;                             % ���ջ�����ϵ�� dB;
S = 350e3*350e3;                    %̽����� (�Զ���)
dS = Rs*lambda/48*Rs*lambda/3;      %ÿ����λ�����

A = 3*48;                           %���߿׾�m
SNR = ( Pa*A^2*10^(RCS/10)*Ti ) / ( 4*pi*lambda^2*k*T0*10^(Fn/10)*Rs^4*10^(Ls/10) );
SNR = 10*log10(SNR);

%% ���ʿ׾���
PA = 4*pi*k*T0*Fn*Rs^2*10^(SNR/10)*S*sin(graze/180*pi)/10^(RCS/10);
PA = 10*log10(PA);

%% ����ͼ 
Az0 = 91.42;                 %% ����������
El0 = 50.18;
Ir1 = taylorwin(Nr_az).';    %%����̩�ռ�Ȩ
Ir2 = taylorwin(Nr_el).';    %%����̩�ռ�Ȩ
It1 = taylorwin(Nr_az,4,-23).';    %%����̩�ռ�Ȩ
It2 = taylorwin(Nr_el,4,-23).';    %%����̩�ռ�Ȩ
% It1 = ones(1,Nr_az);    %%������ȼ�Ȩ
% It2 = ones(1,Nr_el);    %%������ȼ�Ȩ
Fr = zeros(length(Az),length(El));  %���շ���ͼ
Ft = zeros(length(Az),length(El));  %���䷽��ͼ
for i=1:length(Az)
    i
    for j = 1:length(El)% 
        Fr(i,j) = fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Ir1,Ir2);
        Ft(i,j) = fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,It1,It2);
    end
end
Fr = 10*log10(abs(Fr*25));
Fr(Fr<-40)=-40;
Ft = 10*log10(abs(Ft*25));
Ft(Ft<-40)=-40;
% % ���գ���ͼ
% [label_az, label_el] = meshgrid(El,Az);
% figure()
% mesh(label_az,label_el,Fr)%10*log10
% xlabel('������/deg')
% ylabel('��λ��/deg')
% zlabel('����/dB')
% title('���շ���ͼ')
% [maxAz_r,maxEl_r] = find(Fr == max(max(Fr)));
% figure()
% plot(Az,Fr(:,maxEl_r))
% title('���շ�λ����ͼ')
% xlabel('��λ��/deg')
% ylabel('����/dB')
% figure()
% plot(El,Fr(maxAz_r,:))
% title('���ո�������ͼ')
% xlabel('������/deg')
% ylabel('����/dB')
% % % ���䣬��ͼ
% [label_az, label_el] = meshgrid(El,Az);
% figure()
% mesh(label_az,label_el,Ft)%10*log10
% xlabel('������/deg')
% ylabel('��λ��/deg')
% zlabel('����/dB')
% title('���䷽��ͼ')
% figure()
% imagesc(El,Az,Ft)%10*log10
% xlabel('������/deg')
% ylabel('��λ��/deg')
% title('���䷽��ͼ')
% [maxAz_t,maxEl_t] = find(Ft == max(max(Ft)));
% figure()
% plot(Az,Ft(:,maxEl_t))
% title('���䷽λ����ͼ')
% xlabel('��λ��/deg')
% ylabel('����/dB')
% figure()
% plot(El,Ft(maxAz_t,:))
% hold on 
% plot(El,Ft(maxAz_r,:),'r')
% title('���丩������ͼ')
% xlabel('������/deg')
% ylabel('����/dB')
% figure()
% index = find(Ft<10);
% label_az_t = label_az;
% label_az_t(index)=nan;
% label_el_t = label_el;
% label_el_t(index)=nan;
% Ft_t = Ft;
% Ft_t(index)=nan;
% contour(label_az_t,label_el_t,Ft_t)
%% ����ģ��
Rsmax = fun_Rsmax(H);                   %���̽��б��m
r = c*Tp/2/D;                           %����ֱ���m
Rsu=c/(2*fr);                           %�����ģ��б��m
Nk1 = floor((Rsmax-Rs)/Rsu);            %ǰ��ģ����
Nk2 = floor((Rs-H)/Rsu);                %����ģ����
Rsk1 = Rs + (0:Nk1)*Rsu;                %ǰ��ģ��б��m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;             %����ģ��б��m
Rsk = [Rsk2,Rsk1];                      %ģ��б��m
Rk = fun_Rs2R(H,Rsk);                   %ģ���ؾ�m
grazek =  fun_GrazeAngle(H,Rk,Rsk);     %��ģ�����������deg
sk = r./cos(grazek/180*pi)*2*pi.*Rk/Nc; % �Ӳ������m^2
% sak = sk*Nc;                            %�Ӳ������m^2
%% ����ɢ��ϵ��
%<<���Ӳ������л���Ԥ���״����þ������>>
for i = 1:length(grazek)
    sigma0(1,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,1));
    sigma0(2,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,2));
    sigma0(3,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,3));
    sigma0(4,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,4));  
    sigma0(5,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,5,2));
end
% figure()
% hold on
% plot(grazek,sigma0(1,:),'r')
% plot(grazek,sigma0(2,:),'g')
% plot(grazek,sigma0(3,:),'b')
% plot(grazek,sigma0(4,:),'k')
% plot(grazek,sigma0(5,:),'c')
% grid on
% box on
% legend('ɳĮ','ũ��','����','��ɽ','����')
% axis([0,76,-80,20]);
% xlabel('�����/deg')
% ylabel('RCS(\sigma_0)/dB')
% title('�Ӳ�����ɢ��ϵ��������ǵı仯')

