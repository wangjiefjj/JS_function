clc
clear 
close all
%%��������
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
rou = 0.90;  %%Э����������ɵĳ�������
sigma_t =0.1;
%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SCRout=0:1:30; % ���SCR
SCRnum=10.^(SCRout/10);
CNRout = 30; %%�����
CNRnum=10.^(CNRout/10);
cos2=0.9;
PFA=1e-2;% PFA=1e-4;

MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e3;
L=round(n*N); 
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��

%% ����Э����
%%�Ӳ�Э����1
Rc = CNRnum*fun_rho(rou,N,1,0.1);
R = Rc+eye(N);
%%�Ӳ�Э����2
% fc = 0;
% sigmaf = 0.03; %%�Ӳ���չ��
% rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
% Rc = CNRnum * toeplitz(rc);
% R = Rc+eye(N);
%%%%%%%%%
R_half=R^0.5;
iR=inv(R);
tic
% h = waitbar(1,'Please wait...');
parfor i = 1:MonteCarloPfa

    warning('off')
%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); % �����źŽ������Ӳ������� 
    %%%%%%%%%%%%%%%%%%%%%Э�������
        R_SCM = (fun_SCMN(Train));
        R_SCM2 = R_SCM + eye(N);
        R_NSCM = (fun_NSCMN(Train));
        R_MESS = fun_MESS(R_SCM);
        R_SFP = fun_SFP(Train,1);
    %%%%%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tanmf_R(i) = fun_AMF(R,x0,s);
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_AMF(R_SCM,x0,s);
    %%%%%% ANMF_SCM2
    Tanmf_SCM2(i) = fun_AMF(R_SCM2,x0,s);
    %%%%%% ANMF_NSCM
    Tanmf_NSCM(i) = fun_AMF(R_NSCM,x0,s);
    %%%%%% ANMF_MESS
    Tanmf_MESS(i) = fun_AMF(R_MESS,x0,s);
    %%%%%% ANMF_SFP
    Tanmf_SFP(i) = fun_AMF(R_SFP,x0,s);  
end
% close(h)
toc
TANMF_R=sort(Tanmf_R,'descend');
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_SCM2=sort(Tanmf_SCM2,'descend');
TANMF_NSCM=sort(Tanmf_NSCM,'descend');
TANMF_MESS=sort(Tanmf_MESS,'descend');
TANMF_SFP=sort(Tanmf_SFP,'descend');

Th_R = (TANMF_R(floor(MonteCarloPfa*PFA-1))+TANMF_R(floor(MonteCarloPfa*PFA)))/2;
Th_SCM = (TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_SCM2 = (TANMF_SCM2(floor(MonteCarloPfa*PFA-1))+TANMF_SCM2(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM = (TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;
Th_MESS = (TANMF_MESS(floor(MonteCarloPfa*PFA-1))+TANMF_MESS(floor(MonteCarloPfa*PFA)))/2;
Th_SFP = (TANMF_SFP(floor(MonteCarloPfa*PFA-1))+TANMF_SFP(floor(MonteCarloPfa*PFA)))/2;
%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_R=0;
counter_SCM=0;
counter_SCM2=0;
counter_NSCM=0;
counter_MESS=0;
counter_SFP=0;

Pd_R_mc = zeros(1,length(SCRout));
Pd_SCM_mc = zeros(1,length(SCRout));
Pd_SCM2_mc = zeros(1,length(SCRout));
Pd_NSCM_mc = zeros(1,length(SCRout));
Pd_MESS_mc = zeros(1,length(SCRout));
Pd_SFP_mc = zeros(1,length(SCRout));
alpha=sqrt(SCRnum/abs(s'*R*s)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
if str_train=='g'
    alpha=sqrt(SCRnum/abs(s'*Rc*s));
elseif str_train=='p'
    alpha=sqrt(SCRnum*mu/(lambda-1));     
end
% alpha=sqrt(SNRnum*N/abs(s'*s));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SCRout)
    waitbar(m/length(SCRout),h,sprintf([num2str(m/length(SCRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); % �����źŽ������Ӳ�������       
    %%%%%%%%%%%%%%%Э�������
        R_SCM = (fun_SCMN(Train));
        R_SCM2 = R_SCM + eye(N);
        R_NSCM = (fun_NSCMN(Train));
        R_MESS = fun_MESS(R_SCM);
        R_SFP = fun_SFP(Train,1);
        %����ź�
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%OPT
        T_R = fun_AMF(R,x0,s);
        %%%%% ANMF_SCM
        T_SCM = fun_AMF(R_SCM,x0,s);  
        %%%%% ANMF_NSCM
        T_NSCM = fun_AMF(R_NSCM,x0,s); 
        %%%%% ANMF_SCM2
        T_SCM2 = fun_AMF(R_SCM2,x0,s);
        %%%%% ANMF_MESS
        T_MESS = fun_AMF(R_MESS,x0,s);
        %%%%% ANMF_SFP
        T_SFP = fun_AMF(R_SFP,x0,s);  
%         %%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T_R>Th_R;                  counter_R=counter_R+1;          end   
        if T_SCM>Th_SCM;              counter_SCM=counter_SCM+1;      end
        if T_SCM2>Th_SCM2;            counter_SCM2=counter_SCM2+1;    end
        if T_NSCM>Th_NSCM;            counter_NSCM=counter_NSCM+1;    end
        if T_MESS>Th_MESS;            counter_MESS=counter_MESS+1;    end    
        if T_SFP>Th_SFP;              counter_SFP=counter_SFP+1;      end
    end
    Pd_R_mc(m)=counter_R/MonteCarloPd;                counter_R=0;
    Pd_SCM_mc(m)=counter_SCM/MonteCarloPd;            counter_SCM=0;
    Pd_SCM2_mc(m)=counter_SCM2/MonteCarloPd;          counter_SCM2=0;
    Pd_NSCM_mc(m)=counter_NSCM/MonteCarloPd;          counter_NSCM=0;
    Pd_MESS_mc(m)=counter_MESS/MonteCarloPd;          counter_MESS=0;
    Pd_SFP_mc(m)=counter_SFP/MonteCarloPd;            counter_SFP=0;
end
toc
close(h)
figure();
hold on
plot(SCRout,Pd_R_mc,'r.-','linewidth',2)
plot(SCRout,Pd_SCM_mc,'g.-','linewidth',2)
plot(SCRout,Pd_SCM2_mc,'b.-','linewidth',2)
plot(SCRout,Pd_NSCM_mc,'k.-','linewidth',2)
plot(SCRout,Pd_MESS_mc,'c.-','linewidth',2)
plot(SCRout,Pd_SFP_mc,'y.-','linewidth',2)
grid on
box on