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
load Th_4Second_s0.9_PFA3_p.mat
%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fd = -0.5:0.05:0.5;
h = waitbar(0,'Please wait...');
%     alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
if str_train=='g'
   alpha=sqrt(SNRnum/abs(s'*irouR*s));
elseif str_train=='p'
    alpha=sqrt(SNRnum*mu/(lambda-1));     
end
Pd_R_mc = zeros(1,length(fd));
Pd_CC_mc = zeros(1,length(fd));
Pd_E_mc = zeros(1,length(fd));
Pd_ECC_mc = zeros(1,length(fd));
Pd_LogM_mc = zeros(1,length(fd));
Pd_LogCC_mc = zeros(1,length(fd));
Pd_P_mc = zeros(1,length(fd));
Pd_PCC_mc = zeros(1,length(fd));
Pd_SFP_mc = zeros(1,length(fd));
tic
for i_t = 1:length(fd)
    waitbar(i_t/length(fd),h,sprintf([num2str(i_t/length(fd)*100),'%%']));
    s = exp(-1i*2*pi*nn*fd(i_t)).'/sqrt(N); %%%%%% 系统导向矢量
    counter_R=0;
    counter_CC=0;
    counter_E=0;
    counter_ECC=0;
    counter_LogM=0;
    counter_LogCC=0;
    counter_P=0;
    counter_PCC=0;
    counter_SFP=0;

%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        Train = awgn(Train,CNR);
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
        x0 = awgn(x0,CNR);
% % %       RKA
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
        R_KA = zeros(size(rouR));
        R_KA = rouR.*(t*t');        
% %         协方差估计
        R_CC = fun_CC(Train,fun_SCMN(Train),R_KA);
        R_E = fun_RPowerEMean(Train,1,4);
        R_ECC = fun_PowerCC(Train,R_KA,1,4);
        R_LogM = fun_RLogEMean(Train,4);
        R_LogCC = fun_LogCC_new(Train,R_KA,4);
        R_P = fun_RPowerEMean(Train,-1,4);
        R_PCC = fun_PowerCC(Train,R_KA,-1,4);
        R_SFP = fun_SFP(Train,1);
        %检测信号
        x0=alpha*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%OPT
        T_R = fun_ANMF(rouR,x0,s);
%         %%%%% ANMF_CC
        T_CC = fun_ANMF(R_CC,x0,s);
%         %%%%%% ANMF_E
        T_E = fun_ANMF(R_E,x0,s);
%         %%%%% ANMF_NSCM
        T_ECC = fun_ANMF(R_ECC,x0,s); 
% %         %%% ANMF_LogM
        T_LogM = fun_ANMF(R_LogM,x0,s);
        if T_LogM>1
            T_LogM = 1;
        end
%         %%%% ANMF_LogCC
        T_LogCC = fun_ANMF(R_LogCC,x0,s);
        %%%%% ANMF_P
        T_P = fun_ANMF(R_P,x0,s);
%         %%%%% ANMF_PCC
        T_PCC = fun_ANMF(R_PCC,x0,s);
%         %%%%% ANMF_SFP
        T_SFP = fun_ANMF(R_SFP,x0,s);  
%         %%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T_R>Th_R;                counter_R=counter_R+1;    end   
        if T_CC>Th_CC;              counter_CC=counter_CC+1;    end
        if T_E>Th_E;            counter_E=counter_E+1;    end
        if T_ECC>Th_ECC;            counter_ECC=counter_ECC+1;    end
        if T_LogM>Th_LogM;          counter_LogM=counter_LogM+1;    end
        if T_LogCC>Th_LogCC;        counter_LogCC=counter_LogCC+1;    end
        if T_P>Th_P;                counter_P=counter_P+1;    end
        if T_PCC>Th_PCC;            counter_PCC=counter_PCC+1;    end
        if T_SFP>Th_SFP;            counter_SFP=counter_SFP+1;    end
    end
    Pd_R_mc(i_t)=counter_R/MonteCarloPd;              counter_R=0;
    Pd_CC_mc(i_t)=counter_CC/MonteCarloPd;            counter_CC=0;
    Pd_E_mc(i_t)=counter_E/MonteCarloPd;              counter_E=0;
    Pd_ECC_mc(i_t)=counter_ECC/MonteCarloPd;          counter_ECC=0;
    Pd_LogM_mc(i_t)=counter_LogM/MonteCarloPd;        counter_LogM=0;
    Pd_LogCC_mc(i_t)=counter_LogCC/MonteCarloPd;      counter_LogCC=0;
    Pd_P_mc(i_t)=counter_P/MonteCarloPd;              counter_P=0;
    Pd_PCC_mc(i_t)=counter_PCC/MonteCarloPd;          counter_PCC=0;
    Pd_SFP_mc(i_t)=counter_SFP/MonteCarloPd;          counter_SFP=0;
end
toc
close(h)
figure();
hold on
plot(fd,Pd_R_mc,'c.-','linewidth',1)
plot(fd,Pd_CC_mc,'b-*','linewidth',1)
plot(fd,Pd_E_mc,'k-.*','linewidth',1)
plot(fd,Pd_ECC_mc,'k-*','linewidth',1)
plot(fd,Pd_LogM_mc,'r.-','linewidth',1)
plot(fd,Pd_LogCC_mc,'r-*','linewidth',1)
plot(fd,Pd_P_mc,'g.-','linewidth',1)
plot(fd,Pd_PCC_mc,'g-*','linewidth',1)
plot(fd,Pd_SFP_mc,'c-*','linewidth',1)

h_leg = legend('NMF','ANMF with CC',...
    'ANMF with E','ANMF with ECC','ANMF with LogM','ANMF with LogCC',...
    'ANMF with P','ANMF with PCC','SFP');

% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
% grid on
% grid minor
% str = ['PdFd_P_',num2str(L),'Second','_',str_train,'.mat'];
str = ['PdFd_',num2str(L),'Second','_s',num2str(sigma_t),'_',str_train,'.mat'];
% % % str = ['Pd_E',num2str(L),'Second','_',str_train,'.mat'];
save (str); 
