%% ���� ����״��ʱ����Ӧ�Ӳ����Ƽ���,������ �����µ�
clc
clear
close all
% isRu: ���޾���ģ������1����0
% isRotation��������ת����1����0
isRu = 0;
isRotation = 0;
%% ����
Re = 6373e3;                    %Բ����뾶m
C  = 3e8;                 %����m/s 299792458
%% �״�ϵͳ����
fo = 1250e6;                    %��ƵHz
lambda = C/fo;                  %����
Nr_az = 16;%16                     %���շ�λ�������
Nr_el = 4;%25                     %���ո����������
Nt_az = 16;%16                     %���䷽λ�������
Nt_el = 4;%25                     %���丩���������
Ti = 14.2e-3;                   %פ��ʱ��s
fr = 500;                      %�����ظ�����Hz
Np = 16;%floor(fr*Ti);              %���������� 
Pt = 30e3;                      %��ֵ����W
Tp = 20e-6;                     %����s
Pa = Pt*Tp/(1/fr);              %ƽ�����书��
Ls = 6;                         %ϵͳ���dB
B = 1e6;                        %����Hz
D = Tp*B;                       %��ѹ��
Gt = 40;                        %��������dB
Gr = 23;                        %��������dB
d = lambda/2;                   %��Ԫ���m
d_saz = 0.12*25;                %��λ��������m
L = 48;                         %���߳���m
W = 3;                          %���߿��m
%% �״�
% alpha_r = 30;                   %�״�γ��N,deg
% beta_r = 110;                   %�״ﾭ��E,deg
H = 850e3;                      %�߶�m
Vp = fun_Vp(H);                 %ƽ̨�ٶ�
eta = 70;                       %������,deg
alpha1 = 110;                   %�״ﾭ��E,deg
beta1 =  30;                    %�״�γ��N,deg
beta = Vp*2/fr/d;               %%���ϵ��;
%% Ŀ��
graze = 30;                     %���߷��������deg
alpha2 = 120;                   %Ŀ�꾭��E,deg
beta2  = 26;                    %Ŀ��γ��N,deg
RCS = 5;                        %Ŀ��RCS��dBsm����dB*m^2��
% ����graze�õ��ĵؾ࣬�����߷���Ĳ���ָ��ؾ�
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
%% ����ָ�����
All0 = fun_ComputeAzEl(alpha1,beta1,H,eta,alpha2,beta2,0);

[X1,Y1,Z1] = fun_JWH2XYZ(alpha1,beta1, H+Re); %���Ǿ�γ��תXYZ
[X2,Y2,Z2] = fun_JWH2XYZ(alpha2,beta2, Re);   %Ŀ��㾭γ��תXYZ
Rs = sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2);     %Ŀ��б��
R = fun_Rs2R(H,Rs);                           %Ŀ��ؾ�
% %���㲨��ָ��ĸ�����
% graze0 = fun_GrazeAngle(H,R,Rs);
% El0 = fun_ELAngle(H,R);
% % ���㲨��ָ��ķ�λ��
% % �����״�����ϵ����
% mu = asin(sin(beta1/180*pi)/sin(eta/180*pi));
% [Xr1,Yr1,Zr1] = fun_XYZ2Radar(alpha1,beta1,eta,Re+H, X1,Y1,Z1);
% Az0 = abs(acos(Xr1/Rs/sin(El0/180*pi)))/pi*180;
%% ����ͼ (��������ͼɶ����������㻹Ҫ������)
Iraz = taylorwin(Nr_az).';                   %%��λ����̩�ռ�Ȩ
Irel = taylorwin(Nr_el).';                   %%��������̩�ռ�Ȩ
Itaz = taylorwin(Nt_az).'; %,4,-23         %%��λ����̩�ռ�Ȩ
Itel = taylorwin(Nt_el).'; %,4,-23         %%��������̩�ռ�Ȩ
% It1 = ones(1,Nr_az);                        %%������ȼ�Ȩ
% It2 = ones(1,Nr_el);                        %%������ȼ�Ȩ
Az0 = All0.Az;   %91.42                       %% ����������
El0 = All0.El;   %50.18                       %%������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�����(2.27)deg
% graze0 = acos(sin(El0/180*pi)*(Re+H)/Re)/pi*180; 
%(2.26)����
% R = ( pi/2 - graze0/180*pi - asin( Re/(Re+H) * cos(graze0/180*pi) ) ) * Re; 
% Rs = fun_R2Rs(H,R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fr = zeros(length(Az),length(El));          %���շ���ͼ
Ft = zeros(length(Az),length(El));          %���䷽��ͼ
% % ÿһ����һ����λ�ǣ�ÿһ����һ��������
for i=1:length(Az)
%     i 
    for j = 1:length(El)% 
        Fr(j,i) = fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Iraz,Irel);
%         Fr(j,i) = fun_Ax(Nr_az,Nr_el,d,d, lambda, Az(i),El(j),Az0,El0,Ir1,Ir2);
        Ft(j,i) = fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,EL,Itaz,Itel);
%         Ft(j,i) = fun_Ax(Nt_az,Nt_el,d,d, lambda, Az(i),El(j),Az0,El0,It1,It2);
    end
end
Fr = 10*log10(abs(Fr));
Fr(Fr<-40)=-40;
Ft = 10*log10(abs(Ft));
Ft(Ft<-40)=-40;
% %%���գ���ͼ
figure()
mesh(Az,El,Fr)%10*log10
xlabel('��λ��/deg')
ylabel('������/deg')
zlabel('����/dB')
title('���շ���ͼ')
% hold on
% plot(Az0,El0,'r.','MarkerSize',20)
% [maxAz_r,maxEl_r] = find(Fr == max(max(Fr)));
% figure()
% plot(Az,Fr(:,maxEl_r))
% title('���ո�������ͼ')
% xlabel('����/deg')
% ylabel('��λ��/dB')
% figure()
% plot(El,Fr(maxAz_r,:))
% title('���շ�λ����ͼ')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����Ŀ�꾭γ�Ȼ���������������ͼ
% LL = 200;
% alpha2_t = linspace(alpha2-20,alpha2+20,LL);
% beta2_t = linspace(beta2-20,beta2+20,LL);
% Fr = zeros(LL,LL);
% Ft = zeros(LL,LL);
% for i_a = 1:LL
%     i_a
%     for i_b = 1:LL
%         All(i_a,i_b) = fun_ComputeAzEl(alpha1,beta1,H,eta,alpha2_t(i_a),beta2_t(i_b),0);
% %         Fr(i_a,i_b) = fun_F(Nr_az,Nr_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
% %             All0.Az,All0.El,EL,Ir1,Ir2);
%         Fr(i_a,i_b) = fun_Ax(Nr_az,Nr_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
%             All0.Az,All0.El,Ir1,Ir2);
% %         Ft(i_a,i_b) = fun_F(Nt_az,Nt_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
% %             All0.Az,All0.El,EL,It1,It2);
%         Ft(i_a,i_b) = fun_Ax(Nt_az,Nt_el,d,d, lambda, All(i_a,i_b).Az,All(i_a,i_b).El,...
%             All0.Az,All0.El,It1,It2);
%     end
% end
% 
% 
% Fr = 10*log10(abs(Fr));
% % Fr(Fr<-40)=-40;
% Ft = 10*log10(abs(Ft));
% Ft(Ft<-40)=-40;
% 
% figure
% imagesc(beta2_t,alpha2_t,abs(Fr))
% xlabel('γ��')
% ylabel('����')
% hold on 
% plot(beta1,alpha1,'g.','MarkerSize',20) %%�״�λ��
% hold on 
% plot(beta2,alpha2,'r.','MarkerSize',20) %%������ָ�� 
% % axis([beta2-5,beta2+5,alpha2-10,alpha2+5,-20,25])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    Rsk = [Rsk2,Rsk1];%Rs;%                 %ģ��б��m
elseif isRu == 0
    Rsk = Rs;%
end
Rk = fun_Rs2R(H,Rsk);                       %ģ���ؾ�m
elk = fun_ELAngle(H,Rk);                    %��ģ�����븩����deg
grazek =  fun_GrazeAngle(H,Rk,Rsk);         %��ģ�����������deg
% grazek = acos(sin(elk/180*pi)*(Re+H)/Re)/pi*180; 
sk = r./cos(grazek/180*pi)*2*pi.*Rk/Nc;     % �Ӳ���Ԫ���m^2  (2.36)
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
        Ft_t=abs(fun_F(Nt_az,Nt_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,Itaz,Itel)).^2;
%         Ft_t=abs(fun_Ax(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,Itaz,Itel)).^2;
        Gtin(j,i) = Gt*Ft_t; 
        Fr_t=abs(fun_F(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,EL,Iraz,Irel)).^2;
%         Fr_t=abs(fun_Ax(Nr_az,Nr_el,d,d, lambda, Az(i),elk(j),Az0,El0,Iraz,Irel)).^2;
        Grin(j,i) = Gr*Fr_t;
    end
end

% ����ÿ����λ��(�Ӳ���)��CNR
for i = 1:Nc
    for j = 1:length(elk)
         CNR(j,i) = (Pt*Gtin(j,i)*Grin(j,i)*lambda^2*10.^(sigmac(j)/10)/((4*pi)^3*Pn*10^(Ls/10)*(Rsk(j)^4)'));
    end   
end
CNR = sigma2 * CNR; %
% CNR = 10 * log10(CNR);
index = find(CNR==max(CNR));
% CNR(index)=1000000;
% CNR(1:index-1) = 1e-4;
% CNR(index+1:end) = 1e-4;
% figure
% mesh((CNR))
% figure()
% [X,Y] = meshgrid(Az,elk);
% mesh(X,Y,CNR)
% xlabel('��λ��/deg')
% ylabel('������/deg')
% zlabel('CNR/dB')

%% �Ӳ�Э������㣬����ÿ���Ӳ���Ŀռ��ʱ��Ƶ�ʣ�elk��azk��
azk = Az;
cmj1 = sin(elk/180*pi).' * cos(azk/180*pi);                    %����׶��
fspc1 = d/lambda*cmj1;                                       %�ռ�Ƶ��
cmj2 = sin(elk/180*pi).' * sin(azk/180*pi);                    %����׶��
fspc2 = d/lambda*cmj2;                                       %�ռ�Ƶ��
% ��ƫ�ǣ���ƫ����
if isRotation == 1
    CrabA = fun_CrabAngle( beta1,eta, H);%0;%                      %ƫ����/rad
    CrabM = fun_CrabMagnitude( beta1,eta, H); %1;%                  %ƫ������
elseif isRotation == 0
    CrabA = 0;%                  %ƫ����/rad
    CrabM = 1;%                  %ƫ������
end
omegac = beta * d/lambda*2 * CrabM * (sin(elk/180*pi).'*cos(azk/180*pi+CrabA));   %�Ӳ���1��������
% mod(omegac(index),0.5)
Vc = zeros(Nr_az*Np*Nr_el,Nc);
% Vc = zeros(Nr_az*Np,Nc);
Rc = zeros(Nr_az*Np*Nr_el,Nr_az*Np*Nr_el);
% load Ksic.mat
for j = 1:length(elk)
    j
    for k=1:Nc
        % �ռ䵼��ʸ��.
        a(:,k) = exp(1i*2*pi*fspc1(j,k)*(0:Nr_az-1));    % �ռ䵼��ʸ��,��λ��/sqrt(Nr_az)
        b(:,k) = exp(1i*2*pi*omegac(j,k)*(0:Np-1));     % Temporal Steering Vector/sqrt(Np)
        c(:,k) = exp(1i*2*pi*fspc2(j,k)*(0:Nr_el-1));  % �ռ䵼��ʸ��,������/sqrt(Nr_el)
        Vc(:,k) = kron(b(:,k),kron(c(:,k),a(:,k))); 
    end
    
    A = diag(fliplr(CNR(j,:)));
%     A = repmat(CNR(j,:),Nr_az*Np*Nr_el,1);
    Rc = Rc + Vc*A*Vc';
end
 
Rn = sigma2*eye(Nr_az*Np*Nr_el);
% Rn = sigma2*eye(Nr_az*Np);
Rcn = Rc + Rn;
E=abs(eig(Rcn));
E = 10*log10(sort(E,'descend')).';
figure
hold on
plot(E(1:100),'r')
figure(3)
mesh(abs(Rcn))
%% MVD-plot
L2 = 500;
%%��һ��ʱ��Ƶ�� 
omega = linspace(-0.5,0.5,L2);
fsp1 = 0;%d/lambda/2*cmj; 
fsp2 = d/lambda*sin(El0/180*pi);%d/lambda/2*cmj; 
iRcn = inv(Rcn);
for j = 1:length(omega)  
    a = exp(1i*2*pi*fsp1*(0:Nr_az-1)).';    % �ռ䵼��ʸ��./sqrt(Nr_az)
    c = exp(1i*2*pi*fsp2*(0:Nr_el-1)).';%/sqrt(Nr_el)
    b = exp(1i*2*pi*omega(j)*(0:Np-1)).';    % ʱ�䵼��ʸ��/sqrt(Np)
%      v = kron(b,a);                           % ��ʱ����ʸ��.
     v = kron(b,kron(c,a));
     MVD(j) = abs(v'*iRcn*v)^2;
end
figure()
plot(omega,10*log10(abs(MVD)/max(max(abs(MVD)))))%10*log10(abs(MVD)/max(max(abs(MVD))))
xlabel('��һ��ʱ��Ƶ��')
ylabel('SINR/dB')
% axis([-0.5,0.5,-60,-40])
%% MVD����mesh
% L2 = 200;
% fsp = linspace(-0.5,0.5,L2);
% omega = linspace(-0.5,0.5,L2);              %% ��һ��ʱ��Ƶ�� 
% iRcn = inv(Rcn);
% for i = 1:L2
%     i
%     for j = 1:L2
%         a = exp(1i*2*pi*fsp(i)*(0:Nr_az-1)).';    % �ռ䵼��ʸ��.
%         b = exp(1i*2*pi*omega(j)*(0:Np-1)).';    % ʱ�䵼��ʸ��
%         v = kron(b,a);                           % ��ʱ����ʸ��.
%         MVD(i,j) = abs(v'*iRcn*v)^2;
%     end
% end
% 
% figure()
% [X,Y] = meshgrid(omega,fsp);
% mesh(X,Y,-10*log10(abs(MVD)/max(max(abs(MVD)))))
% xlabel('��һ���ռ�Ƶ��')
% ylabel('��һ��ʱ��Ƶ��')
