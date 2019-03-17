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
SNRnum=10.^(SNRout/10);
PFA=1e-2;% PFA=1e-4;
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
Zhh = sig;
before = 100+N-1; %%ȥǰ��֡��Ϊ����Э����
nn = 0:N-1;
% load Th_4Second19980223_170435_IPIX3.mat
fd = -0.01:0.001:0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%���޼���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_fd = 1:length(fd)
    i_fd/length(fd)
    s = exp(-1i*2*pi*nn*fd(i_fd)).'/sqrt(N); %%%%%% ϵͳ����ʸ��
    Tanmf_CC = zeros(1,MonteCarloPfa);
    Tanmf_E = zeros(1,MonteCarloPfa);
    Tanmf_ECC = zeros(1,MonteCarloPfa);
    Tanmf_P = zeros(1,MonteCarloPfa);
    Tanmf_PCC = zeros(1,MonteCarloPfa);
    Tanmf_LogM = zeros(1,MonteCarloPfa);
    Tanmf_LogCC = zeros(1,MonteCarloPfa);
    Tanmf_SFP = zeros(1,MonteCarloPfa);
    parfor i = 1:MonteCarloPfa  
    warning('off')
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
    % % % %     Э�������
        R_CC = fun_CC(Train,fun_SCMN(Train),R_KA2);
        R_E = fun_RPowerEMean(Train,1,3);
        R_ECC = fun_PowerCC(Train,R_KA1,1,10);
        R_LogM = fun_RLogEMean(Train,3);
        R_LogCC = fun_LogCC_new(Train,R_KA1,10);
        R_P = fun_RPowerEMean(Train,-1,3);
        R_PCC = fun_PowerCC(Train,R_KA1,-1,9);
        R_SFP = fun_SFP(Train,1);
    %     %%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ANMF_CC
        Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
    %     %%%%% ANMF_NSCM
        Tanmf_E(i) = fun_ANMF(R_E,x0,s);
        %%%%% ANMF_NSCM
        Tanmf_ECC(i) = fun_ANMF(R_ECC,x0,s); 
         %%% ANMF_LogM
        Tanmf_LogM(i) = fun_ANMF(R_LogM,x0,s);
        if Tanmf_LogM(i)>1
            Tanmf_LogM(i) = 0;
        end
        %%%%% ANMF_LogCC
        Tanmf_LogCC(i) = fun_ANMF(R_LogCC,x0,s);
        %%%%%% ANMF_P
        Tanmf_P(i) = fun_ANMF(R_P,x0,s);
    %     %%%%%% ANMF_PCC
        Tanmf_PCC(i) = fun_ANMF(R_PCC,x0,s);
    %     %%%%%% ANMF_SFP
        Tanmf_SFP(i) = fun_ANMF(R_SFP,x0,s);  
    end
    % close(h)
    TANMF_CC=sort(Tanmf_CC,'descend');
    TANMF_E=sort(Tanmf_E,'descend');
    TANMF_ECC=sort(Tanmf_ECC,'descend');
    TANMF_LogM=sort(Tanmf_LogM,'descend');
    TANMF_LogCC=sort(Tanmf_LogCC,'descend');
    TANMF_P=sort(Tanmf_P,'descend');
    TANMF_PCC=sort(Tanmf_PCC,'descend');
    TANMF_SFP=sort(Tanmf_SFP,'descend');

    Th_CC(i_fd) = (TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;
    Th_E(i_fd) = (TANMF_E(floor(MonteCarloPfa*PFA-1))+TANMF_E(floor(MonteCarloPfa*PFA)))/2;
    Th_ECC(i_fd) = (TANMF_ECC(floor(MonteCarloPfa*PFA-1))+TANMF_ECC(floor(MonteCarloPfa*PFA)))/2;
    Th_LogM(i_fd) = (TANMF_LogM(floor(MonteCarloPfa*PFA-1))+TANMF_LogM(floor(MonteCarloPfa*PFA)))/2;
    Th_LogCC(i_fd) = (TANMF_LogCC(floor(MonteCarloPfa*PFA-1))+TANMF_LogCC(floor(MonteCarloPfa*PFA)))/2;
    Th_P(i_fd) = (TANMF_P(floor(MonteCarloPfa*PFA-1))+TANMF_P(floor(MonteCarloPfa*PFA)))/2;
    Th_PCC(i_fd) = (TANMF_PCC(floor(MonteCarloPfa*PFA-1))+TANMF_PCC(floor(MonteCarloPfa*PFA)))/2;
    Th_SFP(i_fd) = (TANMF_SFP(floor(MonteCarloPfa*PFA-1))+TANMF_SFP(floor(MonteCarloPfa*PFA)))/2;
end
%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fd = -0.5:0.05:0.5;

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
        R_E = fun_RPowerEMean(Train,1,3);
        R_ECC = fun_PowerCC(Train,R_KA1,1,10);
        R_LogM = fun_RLogEMean(Train,3);
        R_LogCC = fun_LogCC_new(Train,R_KA1,10);
        R_P = fun_RPowerEMean(Train,-1,3);
        R_PCC = fun_PowerCC(Train,R_KA1,-1,9);
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
        if T_CC>Th_CC(i_fd);              counter_CC=counter_CC+1;    end
        if T_E>Th_E(i_fd);            counter_E=counter_E+1;    end
        if T_ECC>Th_ECC(i_fd);            counter_ECC=counter_ECC+1;    end
        if T_LogM>Th_LogM(i_fd);          counter_LogM=counter_LogM+1;    end
        if T_LogCC>Th_LogCC(i_fd);        counter_LogCC=counter_LogCC+1;    end
        if T_P>Th_P(i_fd);                counter_P=counter_P+1;    end
        if T_PCC>Th_PCC(i_fd);            counter_PCC=counter_PCC+1;    end
        if T_SFP>Th_SFP(i_fd);            counter_SFP=counter_SFP+1;    end
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
str = ['PdFd_PFA3_',num2str(L),'Second',str_IPIX_t,'IPIX','.mat'];
save (str); 


