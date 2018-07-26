%% ���������ڲ��������״�ز���������
clear 
close all
clc
isRotation = 0;
%% ���ɷ������
%% �״����
N=8;%ÿ����Ԫ��
M=8;%ÿ����Ԫ��
K=8;%������
lambda=0.5;%����
c=3e8;%����
d=0.25;%��Ԫ���
fr=10000;%�����ظ�Ƶ��
% Rmax=para_data(10);%������þ��� �����ò�����������ݳ������㡿
el0=0;%����ָ������
az0=pi/2;%����ָ��λ��
Inn=30;%���������б�ѩ��Ȩ��dB
Imm=30;%���������б�ѩ��Ȩ��dB
Ikk=30;%ʱ�������б�ѩ��Ȩ��dB
% CNR_dB=para_data(16);%�����dB
% CNR=10^(CNR_dB/10);%�����
thetap=0;%����ͺ���֮��н�
senser_error=0;%��Ԫ���
channel_error=0;%ͨ�����
%channel_error=0.05;%ͨ�����
gamma=1;
xinhao=1;
Q=1; %% �ź�����,MIMO�Ļ�Q=8
%% ƽ̨����
H=872e3;%�߶�
V=fun_Vp(H);%�ٶ�
alpha1 = 30;%ƽ̨γ��
eta = 70;%ƽ̨���
if isRotation==1
    crabA = fun_CrabAngle(alpha1,eta, H); %ƫ����rad
    crabM = fun_CrabMagnitude(alpha1,eta, H);%%ƫ������
else
    crabA = 0; %ƫ����rad
    crabM = 1;%%ƫ������ 
end

%% �Ӳ�
Br=0;%�Ӳ����
flag_clutter=4;%�Ӳ����ͣ��ֱ�Ϊ���Ӳ���ɳĮ��ũ����ꡢ��ɽ
flag_sea=5;%����ȼ�1-5
flag_kind=3;%�Ӳ����ȷֲ����ͣ��ֱ�ΪLog-normal��Weibull��K�ֲ�
flag_noise=0;%�Ƿ�������
Rs=eye(M);
% Rs=eye(M)+tril(0.8*rand(M,M),-1)+triu(0.8*rand(M,M),1);
% Rs=eye(M)+tril(0.05*rand(M,M),-1)+triu(0.05*rand(M,M),1);
%====���ò���
N_theta = 180; % ####��λ�����ַ���
N_psi=180;%׶�����ַ���&& û�õ�������
dpsi=pi/180;%׶�ǲ���&& û�õ�������
N_phi=90;%���������ַ���&& û�õ�������
dphi=pi/90;%�����ǲ���&& û�õ�������
dtheta=pi/180;%��λ�ǲ���&& û�õ�������
% % caa=(channel_error*randn(N,1));ca=caa*caa.';%ͨ����������(�м�������)
% % cpp=(channel_error*pi*randn(N,1));cp=cpp*cpp.';%ͨ������λ���(�м���λ���)
ch_err_vector = (1+rand(N,1)*channel_error).*exp(2i*pi*rand(N,1)*channel_error); % ####ͨ���������������Ϊ��1Ϊ��ֵ�ĸ�˹�ֲ�����λΪ0~2�еľ��ȷֲ���
if Inn==0 In=ones(N,1);else In=chebwin(N,Inn);end %�����б�ѩ��Ȩֵ��������
if Imm==0 Im=ones(M,1);else Im=chebwin(M,Imm);end%�����б�ѩ��Ȩֵ��������
if Ikk==0 Ik=ones(K,1);else Ik=chebwin(K,Ikk);end%�����б�ѩ��Ȩֵ��������

%====���㲻ͬ���뻷б��Rl�͸ߵͽ�phil
Re=6378e3;%��������
Rmax=sqrt(H*(2*Re+H));%�״�������þ���
Ru=c/(2*fr);%���ģ������
% ���㣨����ģ�����۵�����L
L=floor((Rmax-H)/Ru)+1;
Rl=H+(1:L)*Ru;%####���ػ��߶Ⱦ�������۵���Ӧ��ʵ�ʾ���    0:L������ʹ�ô�H���봦��ʼ����
phil=asin(H./Rl);%��ͬ�۵��ĳ�ʼ������
R_percell = 150; % һ�����������������ڿռ�ľ���
N_pcell_max = floor((Ru-H)/R_percell); % ��ģ�����뷶Χ���ж��ٸ�������.ע��Ҫ��֤��Ҫ��������ģ����Χ
N_pcell = 30; % ####�����Ÿ�����ע�����ܴ���N_pcell_max
%��ȡ����������G��Ҫ�ļ���
if flag_clutter==1
    F1=4*10^(-7)*10^(0.6*(flag_sea+1));B=pi/2;belt_0=2.44*(flag_sea+1)^1.08/57.29;%��ͬ�����º��Ӳ�ϵ��
elseif  flag_clutter==2
    A=0.00126;B=pi/2;belt_0=0.14;%ɳĮ�����ϵ��
elseif  flag_clutter==3
    A=0.004;B=pi/2;belt_0=0.2;belt_c=1;%ũ�������ϵ��
elseif  flag_clutter==4
    A=0.0126;B=pi/2;belt_0=0.4;belt_c=1;%���������ϵ��
elseif  flag_clutter==5
    A=0.04;B=1.24;belt_0=0.5;belt_c=1;%��ɽ�����ϵ��
end
Pt=180e3;%���书��
Scale=lambda*sqrt(Pt/(4*pi)^3);%G�е�Pt,wavelength,(4*pi)^3
%�Ӳ��������ģ��
if flag_kind==3
    clutter_mean=3;clutter_var=7;%K �ֲ�ʱ��U��V(Ҫ��������)
else
    clutter_mean=0.5;clutter_var=1.2;%log-Normal��weibull�Ӳ���ֵ������
end

%% �Ӳ���Ƶ�׷ֲ�������1Ϊ���շ��ٷ��棻����2Ϊ���ٷֱȷ��档�������ַ�����VmΪ0�򲻲����Ӳ��ڲ��˶�����0ֵ������Ӳ��ڲ��˶���

%��������1:���շ��ٷ���������£������Ƿ��٣�
theta_3dB = 0.88*lambda/((N-1)*d);%����ˮƽ��3dB���(����)
Vm = 0; % ### ����������Ӳ��ڲ��˶�����Ҫ��������Ϊ0�������Ҫ����ICM���������0ֵ����Ϊ����ǿ�ƶ�delt_v��ֵ  %�ء����Ӳ�ʱ�ķ���
% % if flag_clutter==1%��ʽ(2.16)
% %     delta_v1=0.101*Vm;    % ���Ӳ�
% % else
% %     delta_v1=0.0066*Vm;  % ���Ӳ�
% % end
% % delta_v2=V*theta_3dB/(2*sqrt(2*log(2)));%��ʽ(2.17)
% % delt_v=sqrt(delta_v1^2+delta_v2^2);%��ʽ(2.18)
% % % delt_v=sqrt(delta_v1^2);%��ʽ(2.18)
delt_v = 0.5; % ����Ϊ�Ӳ��ڲ��˶��ٶȣ�������������÷�ʽ����ȡֱ��ǿ�Ƹ�ֵ���Ӳ��ٶȱ�׼�� 0.5m/s
delt_f = 2*delt_v/lambda;%��ʽ(2.15)

%��������2:�����Ӳ�����ٷֱȷ������
% % % theta_3dB=0.88*lemda/((N-1)*d);%����ˮƽ��3dB���(����)
% % % delt_f=0.05*2*V/lemda; % ǰ���0.05��ζ���Ӳ��ڲ��˶��ٶ�Ϊ�ɻ��ٶȵ�0.05

f_3db = 2*sqrt(2*log(2))*delt_f;
deltf = 0.6*f_3db; % �ڲ��˶��ٶȶ�Ӧ��Ƶ���Ŷ���Χ
% sk = delt_f*exp(-(deltf*((1:5)-3)).^2/(2*delt_f^2))/sqrt(2*pi*delt_f^2);%��ʽ(2.19)
sk = exp(-(deltf*((1:5)-3)).^2/(2*delt_f^2));%&& �޸��ˣ���J.Ward���棩
Delt = 2*pi*deltf/fr;
sk_sqrt = sqrt(sk);

%% ��������
kk=1.38e-23;%����������(J/KHz)
T0=290;%�¶�
Bn=70e6;%���ջ����� %% ��ȷ����������û��
Fn=3.5;%���ջ�����ϵ��dB
Pn=kk*T0*Bn*Fn;%��������
Rorigin=25000;%��ʼ���þ���
f0=c/lambda/10^9;%����Ƶ��GHz
u=sqrt(f0)/4.7;
ops=0;% �������
Sys_DOF = Q*N*K;   % ϵͳ���ɶ�
tic
%% �Ӳ���������
% for ncell=1:N_pcell%��һ��ѭ������ͬ������
for ncell=1%��һ��ѭ������ͬ������
    ops=ops+1
    Rcell=Rorigin+(ncell-1)*R_percell;%RcellΪ��ͬ�����Ŷ�Ӧ�Ĳ�ͬ����
    clutter_ceho_temp = zeros(Sys_DOF,1);                        % && ���ڴ洢ÿ�������ŵ��Ӳ�����
    clutter_covariance_temp = zeros(Sys_DOF);                  % && ���ڴ洢ÿ�������ŵ��Ӳ�Э�������
    for l=1:L%�ڶ���ѭ������ͬ�������
        l
        Rlc=Rl(l)+Rcell; % RlcΪ��ͬ�����Ų�ͬ���뻷��Ӧ�ľ���           Rl ���Ѿ������ػ��߶� H
        if  Rlc>=H & Rlc<=Rmax
            sin_thetag=H./Rlc-(Rlc.^2-H^2)./(2*Re.*Rlc);%sin(theta_g)
            if sin_thetag>0 & (1-sin_thetag)>(1e-10)
                theta_g=atan(sin_thetag/sqrt(1-sin_thetag^2+(1e-20)));%���ؽ�
                sin_phil=H/Rlc+(Rlc^2-H^2)/(2*Rlc*Re);%��ͬ����ʱ�ĸ����ǵ�����sin(phe)
                phi_l=atan(sin_phil/sqrt(1-sin_phil^2+(1e-20)));%������
                cos_phil=sqrt(1-sin_phil^2);%�����ǵ�����
                if flag_clutter==1
                    theta_c=asin(lambda/(4*pi*(0.025+0.046*(flag_sea^1.72))));
                    if theta_g<=theta_c
                        F3=(theta_g/theta_c)^1.9;
                    else
                        F3=1;
                    end
                    delt0=F1*F3*sin(theta_g)/lambda+u*cot(belt_0)^2*exp(-tan(B-theta_g)^2/tan(belt_0)^2);
                elseif flag_clutter==2
                    theta_c=asin(lambda/(4*pi*(9.3*belt_0^2.2)));
                    if theta_g<=theta_c
                        F3=theta_g/theta_c;
                    else
                        F3=1;
                    end
                    delt0=A*F3*sin(theta_g)/lambda+u*cot(belt_0)^2*exp(-tan(B-theta_g)^2/tan(belt_0)^2);
                else
                    delt0=A*sin(theta_g)/lambda+u*cot(belt_0)^2*exp(-tan(B-theta_g)^2/tan(belt_0)^2);
                end
                sr=Rlc*theta_3dB*R_percell/sqrt(2);
                deltc=delt0*sr;%�Ӳ���Ԫ�״�����
                for ii=1:N_theta%������ѭ������ͬ��λ��
                    if  flag_kind==1
                        al=lognrnd(clutter_mean,clutter_var);
                    elseif  flag_kind==2
                        al=weibell(clutter_var,clutter_mean,0);
                    elseif  flag_kind==3
                        al=-1.9045;%kk_gamm(clutter_mean,clutter_var,0);
                    end
                    if Vm==0
                        X=ones(1,K);
                    else
                        X=cell_f(sk_sqrt,K,Delt);
                    end
                    theta=(ii-1)*(pi/N_theta);
                    ws=2*pi*d*cos(theta)*cos_phil/lambda; % ����Ƶ��
                    wt=4*pi*V*crabM*cos(theta+thetap+crabA)*cos_phil/(lambda*fr); % ʱ��Ƶ��
                    St=exp(1j*wt*(0:K-1)).';
                    Ssr=exp(1j*ws*(0:N-1)).';
                    Sst=exp(1j*ws*gamma*(0:M-1)).';
                    if xinhao==1%%�����  gamma=1
                        Us=[1;zeros(M-1,1)];
                        Ut=ones(M,1);
                    elseif xinhao==2%MIMO,gamma>=1
                        Us=eye(M);
                        Ut=eye(M);
                    elseif xinhao==3
                        Us=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
                        Ut=[0.7 0 0 0;0.7 0 0 0;0 0.7 0 0;0 0.7 0 0;0 0 0.7 0;0 0 0.7 0;0 0 0 0.7;0 0 0 0.7];
                    elseif xinhao==4
                        Us=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
                        Ut=[1 0 0 0;1 0 0 0;0 1 0 0;0 1 0 0;0 0 1 0;0 0 1 0;0 0 0 1;0 0 0 1];
                    end
                    %C=Scale*lognrnd(0,1.2,1,1)/Rlc^2*kron(Us.'*Rs*Us*(Ut.')*(Im.*Sst),kron(St,Ssr));
                    C=Scale*sqrt(deltc)/Rlc^2*kron(Us.'*Rs.'*Us*Ut.'*Sst,kron(St,Ssr));
                    clutter_ceho_temp = clutter_ceho_temp+C;                               % && ���ڴ洢ÿ�������ŵ��Ӳ�����
                    clutter_covariance_temp = clutter_covariance_temp+C*C';         % && ���ڴ洢ÿ�������ŵ��Ӳ�Э�������
                end
            end
        end
    end
    clutter_data(:,ncell) = clutter_ceho_temp; % �������Ӳ����� û������
    clutter_covariance{ncell} = clutter_covariance_temp; % �������Ӳ�Э�������
    clear  clutter_ceho_temp   clutter_covariance_temp;
end
toc
parameter=[fr M N K lambda V H el0 az0 crabA thetap];
save clutterE.mat clutter_data parameter clutter_covariance
ch_num = N; % ��Ԫ��
p_num = K; % ������
inf_str = '����Ԫ������';
% save SL_ULA_8C_8P_300R_IDEAL.mat    clutter_echo    clutter_covariance    inf_str
% save SL_ULA_8C_8P_300R_5.0ICM.mat    clutter_echo    clutter_covariance    inf_str
% save SL_ULA_8C_8P_300R_�ٷ�֮��ICM.mat    clutter_echo    clutter_covariance    inf_str
% save SL_ULA_8C_8P_300R_0.05CE.mat    ch_err_vector    clutter_echo    clutter_covariance    inf_str

%����Э��������ֵD����ͼ
CNR_dB = 60;
CNR = 10^(CNR_dB/10);
R_real = clutter_covariance{1};
%R_real =clutter_echo*clutter_echo'/100;

R = CNR*R_real./sum(eig(R_real)/Sys_DOF);
R = R+1*eye(Sys_DOF);
DD = eig(R);
DD = sort(abs(DD),'descend');
D_db = db(DD)/2;
figure;plot(D_db,'*-');grid on;title('�Ӳ�Э��������������');

%�����Ӳ�������P����ͼ
fd = linspace(-1,1,180);
fs = linspace(-1,1,180);
R_inv=inv(R);
for p=1:length(fs)
    Ssr=exp(1j*pi*(0:N-1)*fs(p)).';
    Sst=exp(1j*pi*gamma*(0:M-1)*fs(p)).';
    for q=1:length(fd)
        St=exp(1j*pi*(0:K-1)*fd(q)).';
        S=kron(Ut.'*Sst,kron(St,Ssr));
        S = S/norm(S);
        P(p,q)=1/abs(S'*R_inv*S);
    end
end
P_db = db(P)/2;
P_db =P_db -max(max(P_db ));
% surf(P_db);shading interp;
figure;
axes('Position',[0.1,0.1,0.8,0.8]);
mesh(fd,fs,P_db);
set(gca,'FontSize',32);
xlabel('2fd/fr','FontSize',32);ylabel('cos\psi_{r}','FontSize',32);zlabel('P/dB','FontSize',32);  %�����Ӳ�������ͼ
axes('Position',[0.6,0.6,0.4,0.4])
contourf(fd,fs,P_db);   %�����Ӳ������׵ȸ���ͼ
set(gca,'xtick',[],'ytick',[]);
axis square


%%
% %% �Ӳ������׷���
%
% %����Э��������ֵD����ͼ
% R =clutter_covariance{1};
% R = CNR*R./sum(eig(R)/Sys_DOF);
% R = R+1*eye(Sys_DOF);
% % R =Clutter*Clutter'/150;
% % CNR_dB = 60;
% % CNR = 10^(CNR_dB/10);
% % %R_real =Clutter_ch*Clutter_ch'/30;
% % %R_real =Clutter*Clutter'/30;
% % R_real = R;
% % R = CNR*R_real./sum(eig(R_real)/Sys_DOF);
% % R = R+1*eye(Sys_DOF);
%
% DD= eig(R);
% DD = sort(abs(DD),'descend');
% D_db = db(DD)/2;
% figure;plot(D_db'-max(D_db),'*-');grid on;title('�Ӳ�Э��������������');
%
% %�����Ӳ�������P����ͼ
% Nfs=180;
% fs = (-Nfs:1:Nfs)/Nfs;
% Nfd=180;
% fd = (-Nfd:1:Nfd)/Nfd;
% %fd=2*cos(phi0)*cos(pi*fs+natenna_angle+crab_theta)*V/lambda/fr;
% R_inv=inv(R);
% for p=1:length(fs)
%     s_ph_azimuth=exp(-j*pi*fs(p)*(0:N-1));                         %��Ƶ��λά��λ
%     s_ph_pitch=exp(-j*pi*fs(p)*(0:0));                                 %��Ƶ����ά��λ
%     Ss=kron(s_ph_pitch,s_ph_azimuth);                                                                   %�ϳɵĿ�Ƶ��λ
%     for q=1:length(fd)
%         St=exp(j*pi*(0:K-1).'*fd(q));
%         S=kron(St,Ss.');
%         S = S/norm(S);
%         P(p,q)=1/abs(S'*R_inv*S);
%     end
% end
% P_db = db(P)/2;
% P_db =P_db -max(max(P_db ));
% %P_db =P_db ;
% figure; %surf(P_db);shading interp;
% axis([-0.5,0.5,-0.5,0.5,-60,0]);
% angle=180*acos(fs)/pi;
%  %mesh(fd,fs,P_db);
% surf(fd,fs,P_db);shading interp;lighting gouraud;colorbar;axis tight;
% title('�Ӳ�������');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB');  %�����Ӳ�������ͼmesh(fd,fs,P_db);cos(\psi)
% figure;contour(fd,fs,P_db,20);   %�����Ӳ������׵ȸ���ͼ
% title('�Ӳ������׵ȸ���');xlabel('2fd/fr');ylabel('cos(\psi)');zlabel('P/dB');
% % P_db = db(P)/2;
% % P_db =P_db -max(max(P_db ));
% % % surf(P_db);shading interp;
% % figure;
% % axes('Position',[0.1,0.1,0.8,0.8]);
% % mesh(fd,fs,P_db);
% % set(gca,'FontSize',32);
% % xlabel('2fd/fr','FontSize',32);ylabel('cos\psi_{r}','FontSize',32);zlabel('P/dB','FontSize',32);  %�����Ӳ�������ͼ
% % axes('Position',[0.6,0.6,0.4,0.4])
% % contourf(fd,fs,P_db);   %�����Ӳ������׵ȸ���ͼ
% % set(gca,'xtick',[],'ytick',[]);
% % axis square
%
% %
% %