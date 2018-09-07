%%%%��������
clc
clear 
close all
%19980223_170435_ANTSTEP.CDF range = 8;
%19980204_224024_ANTSTEP.CDF range = 17;
str_IPIX = '19980223_170435_ANTSTEP.CDF';
str_IPIX_t = str_IPIX(1:16);
[sig,Range,matFile]=fun_Data_process(8,str_IPIX);
load(matFile)
%lambda  %2.4072��19980223_170435��%1.1967(19980204_224024)
%mu      %1.3600��19980223_170435��%1.3180(19980204_224024)
lambda =  2.4072;   
mu = 1.3600;       
%%%%��������
n = 0.5; %����������%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=10; % ���SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPd=1e3;
L=round(n*N); 
Zhh = sig;
before = 80+N-1; %%ȥǰ��֡��Ϊ����Э����

load Th_4Second19980223_170435_IPIX.mat

%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fd = -0.5:0.05:0.5;
nn = 0:N-1;
alpha=sqrt(SNRnum*mu/(lambda-1));  

counter_CC=0;
counter_E=0;
counter_ECC=0;
counter_LogM=0;
counter_LogCC=0;
counter_P=0;
counter_PCC=0;
counter_SFP=0;


Pd_CC_mc = zeros(1,length(fd));
Pd_E_mc = zeros(1,length(fd));
Pd_ECC_mc = zeros(1,length(fd));
Pd_LogM_mc = zeros(1,length(fd));
Pd_LogCC_mc = zeros(1,length(fd));
Pd_P_mc = zeros(1,length(fd));
Pd_PCC_mc = zeros(1,length(fd));
Pd_SFP_mc = zeros(1,length(fd));

      
h = waitbar(0,'Please wait...');
tic
for i_fd=1:length(fd)
    waitbar(i_fd/length(fd),h,sprintf(['�����ʼ���: ', num2str(i_fd/length(fd)*100),'%%']));
    s = exp(-1i*2*pi*nn*fd(i_fd)).'/sqrt(N); %%%%%% ϵͳ����ʸ��
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%0%
        index_t1 = ceil(rand()*(M-10000))+2000;
        R_KA1 = zeros(N,N);
        R_KA2 = zeros(N,N);
        for ii = 1:before-N+1
            x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
            R_KA1 = R_KA1+fun_NSCMN(x_tt)/(before-N+1);
            R_KA2 = R_KA2+x_tt*x_tt'/(before-N+1);
        end
%         R_KA = eye(N);
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
        Train = [Train1,Train2];%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = Zhh(index_t1:index_t1+7,Range) ; % �����źŽ������Ӳ�������
        %%%����ź�
        x0=alpha*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        R_CC = fun_CC(Train,fun_SCMN(Train),R_KA2);
        R_E = fun_RPowerEMean(Train,1,4);
        R_ECC = fun_PowerCC(Train,R_KA1,1,4);
        R_LogM = fun_RLogEMean(Train,4);
        R_LogCC = fun_LogCC_new(Train,R_KA1,4);
        R_P = fun_RPowerEMean(Train,-1,4);
        R_PCC = fun_PowerCC(Train,R_KA1,-1,4);
        R_SFP = fun_SFP(Train,1);


        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%% ANMF_CC
        T_CC = fun_ANMF(R_CC,x0,s);
%         %%%%%% ANMF_NSCM
        T_E = fun_ANMF(R_E,x0,s); 
%         %%%%%% ANMF_NSCM
        T_ECC = fun_ANMF(R_ECC,x0,s); 
%         %%%%%% ANMF_LogM
        T_LogM = fun_ANMF(R_LogM,x0,s);
%         %%%%%% ANMF_LogCC
        T_LogCC = fun_ANMF(R_LogCC,x0,s);
%         %%%%%% ANMF_P
        T_P = fun_ANMF(R_P,x0,s);
        %%%%%% ANMF_PCC
        T_PCC = fun_ANMF(R_PCC,x0,s);
%         %%%%%% ANMF_SFP
        T_SFP = fun_ANMF(R_SFP,x0,s); 
%         %%%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T_CC>Th_CC;              counter_CC=counter_CC+1;    end
        if T_E>Th_E;            counter_E=counter_E+1;    end
        if T_ECC>Th_ECC;            counter_ECC=counter_ECC+1;    end
        if T_LogM>Th_LogM;          counter_LogM=counter_LogM+1;    end
        if T_LogCC>Th_LogCC;        counter_LogCC=counter_LogCC+1;    end
        if T_P>Th_P;                counter_P=counter_P+1;    end
        if T_PCC>Th_PCC;            counter_PCC=counter_PCC+1;    end
        if T_SFP>Th_SFP;            counter_SFP=counter_SFP+1;    end
    end
    Pd_CC_mc(i_fd)=counter_CC/MonteCarloPd;            counter_CC=0;
    Pd_E_mc(i_fd)=counter_E/MonteCarloPd;          counter_E=0;
    Pd_ECC_mc(i_fd)=counter_ECC/MonteCarloPd;          counter_ECC=0;
    Pd_LogM_mc(i_fd)=counter_LogM/MonteCarloPd;        counter_LogM=0;
    Pd_LogCC_mc(i_fd)=counter_LogCC/MonteCarloPd;      counter_LogCC=0;
    Pd_P_mc(i_fd)=counter_P/MonteCarloPd;              counter_P=0;
    Pd_PCC_mc(i_fd)=counter_PCC/MonteCarloPd;          counter_PCC=0;
    Pd_SFP_mc(i_fd)=counter_SFP/MonteCarloPd;          counter_SFP=0;
end
toc
close(h)
figure();
hold on
plot(fd,Pd_CC_mc,'b-*','linewidth',1)
plot(fd,Pd_E_mc,'k-.','linewidth',1)
plot(fd,Pd_ECC_mc,'k-*','linewidth',1)
plot(fd,Pd_LogM_mc,'r.-','linewidth',1)
plot(fd,Pd_LogCC_mc,'r-*','linewidth',1)
plot(fd,Pd_P_mc,'g.-','linewidth',1)
plot(fd,Pd_PCC_mc,'g-*','linewidth',1)
plot(fd,Pd_SFP_mc,'c-*','linewidth',1)

h_leg = legend('ANMF with CC', 'ANMF with E','ANMF with ECC',...
    'ANMF with LogM','ANMF with LogCC','ANMF with P','ANMF with PCC',...
    'ANMF with SFP');

xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
grid on
str = ['PdFd_',num2str(L),'Second',str_IPIX_t,'IPIX','.mat'];
save (str); 


