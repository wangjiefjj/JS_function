%%%������ڲ�ͬ��⻷���µ�����
clc
clear 
close all
%%%%��������
n = 2; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 1;%%%ԽС�Ǹ�˹Խ����
mu = 1;
%%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
%2ÿ����Ԫ����ֵһ����Ϊ���־��Ȼ���
opt_train = 1; 
rou = 0.95;  %%Э����������ɵĳ�������
%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=0:1:35; % ���SNR
SNRnum=10.^(SNRout/10);
CNRout=30;
CNRnum=10.^(CNRout/10); 
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;  
L=round(n*N); 
theta_sig = 0.3;
nn = 0:N-1;
p = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��
%%�Ӳ�Э����
fc=0;
sigmaf = 0.03; %%�Ӳ���չ��
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc = CNRnum * toeplitz(rc);
Rc = Rc+eye(N);
% rouR = fun_rho(rou,N,1,0.1);%%��ʵ���Ӳ�Э����
% irouR = inv(rouR);
tic
% h = waitbar(1,'Please wait...');
%%����
parfor i = 1:MonteCarloPfa
%     warning('off')
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,Rc,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    [x0,tau0] = fun_TrainData(str_train,N,1,Rc,lambda,mu,opt_train); % �����źŽ������Ӳ�������
    %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCM(Train));  
    R_NSCM = (fun_NSCMN(Train));  
    R_AML = fun_AML(Train);
    %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% KGLRT
    Tkglrt(i) = fun_KGLRT(R_SCM,x0,p);
    %%%%% AMF
    Tamf(i) = fun_AMF(R_SCM,x0,p);
    %%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,p);
    %%%%% ANMF_AML
    Tanmf_AML(i) = fun_ANMF(R_AML,x0,p);
end
% close(h)
toc
TKGLRT=sort(Tkglrt,'descend');
TAMF=sort(Tamf,'descend');
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_AML=sort(Tanmf_AML,'descend');

Th_KGLRT = (TKGLRT(floor(MonteCarloPfa*PFA-1))+TKGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_AMF = (TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM=(TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_AML=(TANMF_AML(floor(MonteCarloPfa*PFA-1))+TANMF_AML(floor(MonteCarloPfa*PFA)))/2;

%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_kglrt=0;
counter_amf=0;
counter_nscm=0;
counter_aml=0;

Pd_KGLRT_mc = zeros(1,length(SNRout));
Pd_AMF_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_AML_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,Rc,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        [x0,tau0] = fun_TrainData(str_train,N,1,Rc,lambda,mu,opt_train); % �����źŽ������Ӳ�������
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCM(Train));  
        R_NSCM = (fun_NSCMN(Train));  
        R_AML = fun_AML(Train);      
        %%%����ź�
        x0=alpha(m)*p+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% KGLRT
        Th_kglrt = fun_KGLRT(R_SCM,x0,p);
        %%%%% AMF
        Th_amf = fun_AMF(R_SCM,x0,p);
        %%%%% ANMF_SCM
        Th_scm = fun_ANMF(R_SCM,x0,p);
        %%%%% ANMF_AML
        Thaml = fun_ANMF(R_AML,x0,p);
        %%%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        if Th_kglrt>Th_KGLRT;          counter_kglrt=counter_kglrt+1;        end      
        if Th_amf>Th_AMF;          counter_amf=counter_amf+1;        end                
        if Th_scm>Th_NSCM;       counter_nscm=counter_nscm+1;    end
        if Thaml>Th_AML;       counter_aml=counter_aml+1;    end
    end   
    Pd_KGLRT_mc(m)=counter_kglrt/MonteCarloPd;           counter_kglrt=0;
    Pd_AMF_mc(m)=counter_amf/MonteCarloPd;           counter_amf=0;
    Pd_NSCM_mc(m)=counter_nscm/MonteCarloPd;        counter_nscm=0;
    Pd_AML_mc(m)=counter_aml/MonteCarloPd;        counter_aml=0;
end
toc
close(h)
figure(1);
hold on
plot(SNRout,Pd_KGLRT_mc,'r','linewidth',2)
plot(SNRout,Pd_AMF_mc,'g','linewidth',2)
plot(SNRout,Pd_NSCM_mc,'K-*','linewidth',2)
plot(SNRout,Pd_AML_mc,'k','linewidth',2)
h_leg = legend('KGLRT','AMF','ANMF-SCM','ANMF-AML');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
% str = [str_train,'_',num2str(n),'N','_s',num2str(sigma_t^2),'.mat'];
% save (str); 

