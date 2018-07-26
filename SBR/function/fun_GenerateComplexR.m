function [ Rcn,El0,d,lambda,fr] = fun_GenerateComplexR( isRu,isRotation )
%   ������Ҫ���Ӳ�Э���ͨ����ϸ��������
% isRu: ���޾���ģ������1����0
% isRotation��������ת����1����0
%% ����
Re = 6373e3;                    %Բ����뾶m
C  = 3e8;                       %����m/s 299792458
%% �״�ϵͳ����
fo = 1.25e9;                    %��ƵHz
lambda = C/fo;                  %����
Nr_az = 8;%16                     %���շ�λ�������
Nr_el = 8;%25                     %���ո����������
Nt_az = 8;%16                     %���䷽λ�������
Nt_el = 8;%25                     %���丩���������
Ti = 14.2e-3;                   %פ��ʱ��s
fr = 10000;                      %�����ظ�����Hz
Np = 8;%floor(fr*Ti);              %���������� 
Pt = 30e3;                      %��ֵ����W
Tp = 20e-6;                     %����s
Pa = Pt*Tp/(1/fr);              %ƽ�����书��
Ls = 6;                         %ϵͳ���dB
B = 1e6;                        %����Hz
D = Tp*B;                       %��ѹ��
Gt = 23;                        %��������dB
Gr = 23;                        %��������dB
d = lambda/2;                   %��Ԫ���m(������)
% d_saz = d*25;                %��λ��������m
L = 48;                         %���߳���m
W = 3;                          %���߿��m
%% �״�
% alpha_r = 30;                   %�״�γ��N,deg
% beta_r = 110;                   %�״ﾭ��E,deg
H = 850e3;                        %�߶�m
Vp = fun_Vp(H);                   %ƽ̨�ٶ�
eta = 70;                         %������,deg
alpha1 = 110;                     %�״ﾭ��E,deg
beta1 =  30;                      %�״�γ��N,deg
beta = Vp*2/fr/d;                 %%���ϵ��;
%% Ŀ��
graze = 30;                     %�����deg
alpha2 = 120;                   %Ŀ�꾭��E,deg
beta2  = 25;                    %Ŀ��γ��N,deg
RCS = 5;                        %Ŀ��RCS��dBsm����dB*m^2��
% ����graze�õ��ĵؾ࣬�����߷���Ĳ���ָ��ֱ��Ŀ��ʱ�ĵؾ�
R_graze = ( pi/2 - graze/180*pi - asin( Re/(Re+H) * cos(graze/180*pi) ) ) * Re; %(2.26)����
Rs_graze = fun_R2Rs(H,R_graze);
% ��λ��
Nc = 180;                           %�Ӳ������ 
Az = linspace(0,179,Nc);            %��λ�Ƕ�
El = linspace(0,179,Nc);            %������Ƕ�
% Az = Az/180*pi;
% El = El/180*pi;
EL =  fun_ELAngle(H,R_graze);       %��б�ǣ����߷��߸����ǣ�
LAz = length(Az);
AFt = ones(1,LAz);                  % ��������������
AFr = ones(1,LAz);                  % ��������������

%% �趨����ͼ (��������ͼɶ����������㻹Ҫ������)
Az0 = 91.671396591855710;   %91.42                       %% ����������
El0 = 49.836208625040160;   %50.18                       %%������
%%�����(2.27)deg
graze0 = acos(sin(El0/180*pi)*(Re+H)/Re)/pi*180; 
%(2.26)����
R = ( pi/2 - graze0/180*pi - asin( Re/(Re+H) * cos(graze0/180*pi) ) ) * Re; 
Rs = fun_R2Rs(H,R);
% ����ͼ����
% Ir1 = taylorwin(Nr_az).';                   %%��λ����̩�ռ�Ȩ
% Ir2 = taylorwin(Nr_el).';                   %%��������̩�ռ�Ȩ
Ir1 = ones(1,Nr_az);                      %%���վ��ȼ�Ȩ
Ir2 = ones(1,Nr_el);                      %%���վ��ȼ�Ȩ
% It1 = taylorwin(Nr_az,4,-23).'; %,4,-23           %%��λ����̩�ռ�Ȩ
% It2 = taylorwin(Nr_el,4,-23).'; %,4,-23           %%��������̩�ռ�Ȩ
It1 = ones(1,Nr_az);                      %%������ȼ�Ȩ
It2 = ones(1,Nr_el);                      %%������ȼ�Ȩ
Fr = zeros(length(Az),length(El));          %���շ���ͼ
Ft = zeros(length(Az),length(El));          %���䷽��ͼ
% ÿһ����һ����λ�ǣ�ÿһ����һ��������
% for i=1:length(Az)
% %     i
%     for j = 1:length(El)% 
%         Fr(j,i) = fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Ir1,Ir2);
%         Ft(j,i) = fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,It1,It2);
%     end
% end
% Fr = 10*log10(abs(Fr*25));
% Fr(Fr<-40)=-40;
% Ft = 10*log10(abs(Ft*25));
% Ft(Ft<-40)=-40;
% % %%���գ���ͼ
% % [label_az, label_el] = meshgrid(El,Az);
% % figure()
% % mesh(label_az,label_el,Fr)%10*log10
% % xlabel('��λ��/deg')
% % ylabel('������/deg')
% % zlabel('����/dB')
% % title('���շ���ͼ')
% [maxAz_r,maxEl_r] = find(Fr == max(max(Fr)));
% figure()
% plot(Az,Fr(:,maxEl_r))
% title('���շ�λ����ͼ')
% xlabel('����/deg')
% ylabel('��λ��/dB')
% figure()
% plot(El,Fr(maxAz_r,:))
% title('���ո�������ͼ')
% xlabel('����/dB')
% ylabel('������/deg')
% % % ���䣬��ͼ
% [label_az, label_el] = meshgrid(El,Az);
% figure()
% mesh(label_az,label_el,Ft)%10*log10
% xlabel('��λ��/deg')
% ylabel('������/deg')
% zlabel('����/dB')
% title('���䷽��ͼ')
% figure()
% imagesc(El,Az,Ft)%10*log10
% xlabel('��λ��/deg')
% ylabel('������/deg')
% title('���䷽��ͼ')
% [maxEl_t,maxAz_t] = find(Ft == max(max(Ft)));
% figure()
% plot(Az,Ft(maxEl_t,:))
% title('���䷽λ����ͼ')
% xlabel('��λ��/deg')
% ylabel('����/dB')
% figure()
% plot(El,Ft(:,maxAz_t))
% hold on 
% plot(El,Ft(:,maxAz_r),'r')
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
%% �ն�SNR����
k = 1.3806488e-23;                  % ������������ J/K.
T0 = 300;                           % ��׼���������� Kelvin.
Fn = 3;                             % ���ջ�����ϵ�� dB;
Te = T0*(10^(Fn/10));               % Effective Receiver Temperature in Kelvin.
Ts = 10^(Ls/10)*Te;                 % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                          % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                          % Receiver Noise Power in Watts
sigma2 = 1;                         % Normalized Noise Power in Watts.
S = 350e3*350e3;                    %̽����� (�Զ���)
dS = Rs*lambda/L*Rs*lambda/W;       %ÿ����λ����� ����P17
A = L*W;                           %���߿׾�m
SNR = ( Pa*A^2*10^(RCS/10)*Ti ) / ( 4*pi*lambda^2*k*T0*10^(Fn/10)*Rs^4*10^(Ls/10));%(2.21)
SNR = 10*log10(SNR);

%% ���ʿ׾���
PA = 4*pi*k*T0*Fn*Rs^2*10^(SNR/10)*S*sin(graze/180*pi)/10^(RCS/10);
PA = 10*log10(PA);

%% ����ģ��
Rsmax = fun_Rsmax(H);                       %���̽��б��m
r = C*Tp/2/D;                               %����ֱ���m
Rsu=C/(2*fr);                               %�����ģ��б��m
Nk1 = floor((Rsmax-Rs)/Rsu);                %ǰ��ģ����
Nk2 = floor((Rs-H)/Rsu);                    %����ģ����
Rsk1 = Rs + (0:Nk1)*Rsu;                    %ǰ��ģ��б��m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %����ģ��б��m
if isRu==1
    disp('isRu')
    Rsk = [Rsk2,Rsk1];%Rs;%                 %ģ��б��m
elseif isRu == 0
    disp('noRu')
    Rsk = Rs;%
end
Rk = fun_Rs2R(H,Rsk);                       %ģ���ؾ�m
elk = fun_ELAngle(H,Rk);                    %��ģ�����븩����deg
grazek =  fun_GrazeAngle(H,Rk,Rsk);         %��ģ�����������deg
% grazek = acos(sin(elk/180*pi)*(Re+H)/Re)/pi*180; 
sk = r./cos(grazek/180*pi)*2*pi.*Rk/Nc*10;     % �Ӳ���Ԫ���m^2  (2.36)
% sak = sk*Nc;                              %�Ӳ������m^2
%% ����ɢ��ϵ��
%<<���Ӳ������л���Ԥ���״����þ������>>
for i = 1:length(grazek)
    sigma0(1,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,1));     %ɳĮ
    sigma0(2,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,2));     %ũ��
    sigma0(3,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,3));     %����
    sigma0(4,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,4));     %��ɽ
    sigma0(5,i) = 10*log10(fun_Morchin(grazek(i)/180*pi,fo,5,2));   %����
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
%% ���Ӳ����״��Чɢ����� 
sigmac = 10.^(sigma0(1,:)./10).*sk;                     %�״��Чɢ�����(2.35)
sigmac = 10*log10(sigmac);                              %dB
%% ÿ���Ӳ���Ԫ��CNR
%  ������պͷ�������
for i = 1:Nc
    for j = 1:length(elk)
        Ft_t=abs(fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,It1,It2)).^2;
        Gtin(j,i) = Gt*Ft_t; 
        Fr_t=abs(fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,Ir1,Ir2)).^2;
        Grin(j,i) = Gr*Fr_t;
    end
end
% ����ÿ����λ��(�Ӳ���)��CNR
for i = 1:Nc
    for j = 1:length(elk)
%          CNR(j,i) = sum(Pt.*Gtin(j,i).*Grin(j,i)*lambda^2.*10.^(sigmac.'/10)./((4*pi)^3*Pn*10^(Ls/10).*(Rsk(j).^4).'));
          CNR(j,i) = (Pt*Gtin(j,i)*Grin(j,i)*lambda^2*10.^(sigmac(j)/10)/((4*pi)^3*Pn*10^(Ls/10)*(Rsk(j)^4)'));
    end   
end
CNR = sigma2 * CNR; %
CNR = 10 * log10(CNR);

% figure()
% [X,Y] = meshgrid(Az,elk);
% mesh(X,Y,CNR)
% xlabel('��λ��/deg')
% ylabel('������/deg')
% zlabel('CNR/dB')

%% �Ӳ�Э������㣬����ÿ���Ӳ���Ŀռ��ʱ��Ƶ�ʣ�elk��azk��
azk = Az;
% ��ƫ�ǣ���ƫ����
if isRotation == 1
    disp('isRotation')
    CrabA = fun_CrabAngle( beta1,eta, H);%0;%                      %ƫ����/rad
    CrabM = fun_CrabMagnitude( beta1,eta, H); %1;%                  %ƫ������
elseif isRotation == 0
    disp('noRotation')
    CrabA = 0;%                  %ƫ����/rad
    CrabM = 1;%                  %ƫ������
end


Vc = zeros(Nr_az*Np*Nr_el,Nc);
% Vc = zeros(Nr_az*Np,Nc);
Rc = zeros(Nr_az*Np*Nr_el,Nr_az*Np*Nr_el);
for j = 1:length(elk)
    j
    for k=1:Nc
        cmj = cos(elk(j)/180*pi).' * cos(azk(k)/180*pi);                    %����׶��
        fspc = d/lambda*cmj;                                       %�ռ�Ƶ��
        omegac = beta * d/lambda*2 * CrabM * (cos(elk(j)/180*pi).*cos(azk(k)/180*pi+CrabA));   %�Ӳ���1��������
        % �ռ䵼��ʸ��.
        a = exp(1i*2*pi*fspc*(0:Nr_az-1)).';    % �ռ䵼��ʸ��,��λ��
        c = exp(1i*2*pi*fspc*(0:Nr_el-1)).';    % �ռ䵼��ʸ��,������
        b = exp(1i*2*pi*omegac*(0:Np-1)).';     % Temporal Steering Vector 
        Vc(:,k) = kron(b,kron(c,a));
    end
    A = diag((CNR(j,:)));
    Rc = Rc+Vc*A*Vc';
end
% Rc = Nr_az*Np*Nr_el*Rc/sum(eig(Rc));
Rn = sigma2*eye(Nr_az*Np*Nr_el);
% Rn = sigma2*eye(Nr_az*Np);

Rcn = Rc + Rn;


end

