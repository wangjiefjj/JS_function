%%%%参数设置
clc
clear 
close all
%19980223_170435_ANTSTEP.CDF range = 8;
%19980204_224024_ANTSTEP.CDF range = 17;
str_IPIX = '19980223_170435_ANTSTEP.CDF';
str_IPIX_t = str_IPIX(1:16);
[sig,Range,matFile]=fun_Data_process(8,str_IPIX);
load(matFile)
%lambda  %2.4072（19980223_170435）%1.1967(19980204_224024)
%mu      %1.3600（19980223_170435）%1.3180(19980204_224024)
lambda =  2.4072;   
mu = 1.3600;       
%%%%参数设置
n = 0.5; %几倍的样本%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=-5:1:25; % 输出SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
Zhh = sig;
before = 80+N-1; %%去前几帧作为先验协方差
tic
% h = waitbar(0,'Please wait...');
parfor i = 1:MonteCarloPfa
    warning('off')
%     warning('query')
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
    index_t1 = ceil(rand()*(M-10000))+2000;
    %%R_KA
%     R_KA = zeros(N,N);
%     for ii = 1:before-N+1
%         x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
%         R_KA = R_KA+x_tt*x_tt'/(before-N+1);
%     end
    R_KA = eye(N);
    Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
    Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
    Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
    %%协方差估计
    R_SCM = fun_SCMN(Train);
    R_CC = fun_CC(Train,R_SCM,R_KA);
    R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
    R_NSCM = fun_NSCMN(Train);
    R_ECC = fun_PowerCC(Train,R_KA,1,4);
    R_LogM = fun_RLogEMean(Train,4);
    R_LogCC = fun_LogCC_new(Train,R_KA,4);
    R_P = fun_RPowerEMean(Train,-1,4);
    R_PCC = fun_PowerCC(Train,R_KA,-1,4);
    R_SFP = fun_SFP(Train,1);
    
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,s);
    %%%%%% ANMF_CC
    Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
    %%%%%% ANMF_ML
    Tanmf_ML(i) = fun_ANMF(R_ML,x0,s);
    %%%%%% ANMF_NSCM
    Tanmf_NSCM(i) = fun_ANMF(R_NSCM,x0,s);
    %%%%%% ANMF_NSCM
    Tanmf_ECC(i) = fun_ANMF(R_ECC,x0,s); 
    %%%%%% ANMF_LogM
    Tanmf_LogM(i) = fun_ANMF(R_LogM,x0,s);
    %%%%%% ANMF_LogCC
    Tanmf_LogCC(i) = fun_ANMF(R_LogCC,x0,s);
    %%%%%% ANMF_P
    Tanmf_P(i) = fun_ANMF(R_P,x0,s);
    %%%%%% ANMF_PCC
    Tanmf_PCC(i) = fun_ANMF(R_PCC,x0,s);
    %%%%% ANMF_SFP
    Tanmf_SFP(i) = fun_ANMF(R_SFP,x0,s); 
end
toc

TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_CC=sort(Tanmf_CC,'descend');
TANMF_ML=sort(Tanmf_ML,'descend');
TANMF_NSCM=sort(Tanmf_NSCM,'descend');
TANMF_ECC=sort(Tanmf_ECC,'descend');
TANMF_LogM=sort(Tanmf_LogM,'descend');
TANMF_LogCC=sort(Tanmf_LogCC,'descend');
TANMF_P=sort(Tanmf_P,'descend');
TANMF_PCC=sort(Tanmf_PCC,'descend');
TANMF_SFP=sort(Tanmf_SFP,'descend');


Th_SCM = (TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_CC = (TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;
Th_ML = (TANMF_ML(floor(MonteCarloPfa*PFA-1))+TANMF_ML(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM = (TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;
Th_ECC = (TANMF_ECC(floor(MonteCarloPfa*PFA-1))+TANMF_ECC(floor(MonteCarloPfa*PFA)))/2;
Th_LogM = (TANMF_LogM(floor(MonteCarloPfa*PFA-1))+TANMF_LogM(floor(MonteCarloPfa*PFA)))/2;
Th_LogCC = (TANMF_LogCC(floor(MonteCarloPfa*PFA-1))+TANMF_LogCC(floor(MonteCarloPfa*PFA)))/2;
Th_P = (TANMF_P(floor(MonteCarloPfa*PFA-1))+TANMF_P(floor(MonteCarloPfa*PFA)))/2;
Th_PCC = (TANMF_PCC(floor(MonteCarloPfa*PFA-1))+TANMF_PCC(floor(MonteCarloPfa*PFA)))/2;
Th_SFP = (TANMF_SFP(floor(MonteCarloPfa*PFA-1))+TANMF_SFP(floor(MonteCarloPfa*PFA)))/2;

%%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_R=0;
counter_SCM=0;
counter_CC=0;
counter_ML=0;
counter_NSCM=0;
counter_ECC=0;
counter_LogM=0;
counter_LogCC=0;
counter_P=0;
counter_PCC=0;
counter_SFP=0;


Pd_SCM_mc = zeros(1,length(SNRout));
Pd_CC_mc = zeros(1,length(SNRout));
Pd_ML_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_ECC_mc = zeros(1,length(SNRout));
Pd_LogM_mc = zeros(1,length(SNRout));
Pd_LogCC_mc = zeros(1,length(SNRout));
Pd_P_mc = zeros(1,length(SNRout));
Pd_PCC_mc = zeros(1,length(SNRout));
Pd_SFP_mc = zeros(1,length(SNRout));

alpha=sqrt(SNRnum*mu/(lambda-1));        
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%0%
        index_t1 = ceil(rand()*(M-10000))+2000;
%         %%R_KA
%         R_KA = zeros(N,N);
%         for ii = 1:before-N+1
%             x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
%             R_KA = R_KA+x_tt*x_tt'/(before-N+1);
%         end
        R_KA = eye(N);
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
        Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
        %%%检测信号
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        R_SCM = fun_SCMN(Train);
        R_CC = fun_CC(Train,R_SCM,R_KA);
        R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
        R_NSCM = fun_NSCMN(Train);
        R_ECC = fun_PowerCC(Train,R_KA,1,4);
        R_LogM = fun_RLogEMean(Train,4);
        R_LogCC = fun_LogCC_new(Train,R_KA,4);
        R_P = fun_RPowerEMean(Train,-1,4);
        R_PCC = fun_PowerCC(Train,R_KA,-1,4);
        R_SFP = fun_SFP(Train,1);


        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% ANMF_SCM
        T_SCM = fun_ANMF(R_SCM,x0,s);
        %%%%%% ANMF_CC
        T_CC = fun_ANMF(R_CC,x0,s);
        %%%%%% ANMF_ML
        T_ML = fun_ANMF(R_ML,x0,s);
        %%%%%% ANMF_NSCM
        T_NSCM = fun_ANMF(R_NSCM,x0,s);
        %%%%%% ANMF_NSCM
        T_ECC = fun_ANMF(R_ECC,x0,s); 
        %%%%%% ANMF_LogM
        T_LogM = fun_ANMF(R_LogM,x0,s);
        %%%%%% ANMF_LogCC
        T_LogCC = fun_ANMF(R_LogCC,x0,s);
        %%%%%% ANMF_P
        T_P = fun_ANMF(R_P,x0,s);
        %%%%%% ANMF_PCC
        T_PCC = fun_ANMF(R_PCC,x0,s);
%         %%%%%% ANMF_SFP
        T_SFP = fun_ANMF(R_SFP,x0,s); 
        %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T_SCM>Th_SCM;            counter_SCM=counter_SCM+1;    end   
        if T_CC>Th_CC;              counter_CC=counter_CC+1;    end
        if T_ML>Th_ML;              counter_ML=counter_ML+1;    end
        if T_NSCM>Th_NSCM;          counter_NSCM=counter_NSCM+1;    end   
        if T_ECC>Th_ECC;            counter_ECC=counter_ECC+1;    end
        if T_LogM>Th_LogM;          counter_LogM=counter_LogM+1;    end
        if T_LogCC>Th_LogCC;        counter_LogCC=counter_LogCC+1;    end
        if T_P>Th_P;                counter_P=counter_P+1;    end
        if T_PCC>Th_PCC;            counter_PCC=counter_PCC+1;    end
        if T_SFP>Th_SFP;            counter_SFP=counter_SFP+1;    end
    end
    Pd_SCM_mc(m)=counter_SCM/MonteCarloPd;          counter_SCM=0;
    Pd_CC_mc(m)=counter_CC/MonteCarloPd;            counter_CC=0;
    Pd_ML_mc(m)=counter_ML/MonteCarloPd;            counter_ML=0;
    Pd_NSCM_mc(m)=counter_NSCM/MonteCarloPd;        counter_NSCM=0;
    Pd_ECC_mc(m)=counter_ECC/MonteCarloPd;          counter_ECC=0;
    Pd_LogM_mc(m)=counter_LogM/MonteCarloPd;        counter_LogM=0;
    Pd_LogCC_mc(m)=counter_LogCC/MonteCarloPd;      counter_LogCC=0;
    Pd_P_mc(m)=counter_P/MonteCarloPd;              counter_P=0;
    Pd_PCC_mc(m)=counter_PCC/MonteCarloPd;          counter_PCC=0;
    Pd_SFP_mc(m)=counter_SFP/MonteCarloPd;          counter_SFP=0;
end
toc
close(h)
figure();
hold on
plot(SNRout,Pd_SCM_mc,'b.-','linewidth',1)
plot(SNRout,Pd_CC_mc,'b-*','linewidth',1)
plot(SNRout,Pd_ML_mc,'b->','linewidth',1)
plot(SNRout,Pd_NSCM_mc,'k.-','linewidth',1)
plot(SNRout,Pd_ECC_mc,'k-*','linewidth',1)
plot(SNRout,Pd_LogM_mc,'r.-','linewidth',1)
plot(SNRout,Pd_LogCC_mc,'r-*','linewidth',1)
plot(SNRout,Pd_P_mc,'g.-','linewidth',1)
plot(SNRout,Pd_PCC_mc,'g-*','linewidth',1)
plot(SNRout,Pd_SFP_mc,'c-*','linewidth',1)

h_leg = legend('ANMF with SCM','ANMF with CC','ANMF with ML','ANMF with NSCM',...
    'ANMF with ECC','ANMF with LogM','ANMF with LogCC','ANMF with P','ANMF with PCC','SFP');

xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
grid on
% str = ['Pd_',num2str(L),'Second',str_IPIX_t,'IPIX','.mat'];
% save (str); 


