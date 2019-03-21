%%%GLC-GLRT的ROC曲线
clc
clear 
close all
str_IPIX = '19980223_170435_ANTSTEP.CDF';
str_IPIX_t = str_IPIX(1:16);
[sig,Range,matFile]=fun_Data_process(8,str_IPIX);
load(matFile)
%lambda  %2.4072（19980223_170435）%1.1967(19980204_224024)
%mu      %1.3600（19980223_170435）%1.3180(19980204_224024)
lambda =  1.5383;   
mu = 1.7745;      
Range = 25;        
% %%%%参数设置
n = 1; %几倍的样本
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=[-5,5,10]; % 输出SNR
cos2=0.9;
% PFA=[1e-4,1e-3:1e-2:1e-1+1e-2];
PFA=[1e-2:1e-2:1e-1];
SNRnum=10.^(SNRout/10);
MonteCarloPfa=round(1./PFA*100);
L_Pfa = length(MonteCarloPfa);
MonteCarloPd=1e4;
rou = 0.90;  %%协方差矩阵生成的迟滞因子
rouR = zeros(N,N);  %%真实的杂波协方差
L=round(n*N); 
Zhh = sig;
before = 100+N-1; %%去前几帧作为先验协方差
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% 系统导向矢量
h = waitbar(0,'Please wait...');

for i_Pfa = 1:L_Pfa
    waitbar((i_Pfa/L_Pfa),h,sprintf(['检测门限计算: ',num2str(i_Pfa/L_Pfa*100),'%%']));
    Tanmf_R = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_CC = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_E = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_ECC = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_LogM = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_LogCC = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_P = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_PCC = zeros(MonteCarloPfa(i_Pfa),1);
    Tanmf_SFP = zeros(MonteCarloPfa(i_Pfa),1);
    offset = 60000-MonteCarloPfa(i_Pfa)-N+1;
    parfor i = 1:MonteCarloPfa(i_Pfa)
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
        %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%0%
        index_t1 = i+offset;
        R_KA1 = zeros(N,N);
        R_KA2 = zeros(N,N);
        for ii = 1:before-N+1
            x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
            R_KA1 = R_KA1+fun_NSCMN(x_tt)/(before-N+1);
            R_KA2 = R_KA2+x_tt*x_tt'/(before-N+1);
        end
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
        Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
        R_CC = fun_CC(Train,fun_SCMN(Train),R_KA2);
        R_E = fun_RPowerEMean(Train,1,3);
        R_ECC = fun_PowerCC(Train,R_KA1,1,10);
        R_LogM = fun_RLogEMean(Train,3);
        R_LogCC = fun_LogCC_new(Train,R_KA1,10);
        R_P = fun_RPowerEMean(Train,-1,3);
        R_PCC = fun_PowerCC(Train,R_KA1,-1,9);
        R_SFP = fun_SFP(Train,1);
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ANMF_CC
        Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
        %%%%% ANMF_CC
        Tanmf_E(i) = fun_ANMF(R_E,x0,s);
        %%%%% ANMF_NSCM
        Tanmf_ECC(i) = fun_ANMF(R_ECC,x0,s); 
        %%%%%% ANMF_LogCC
        Tanmf_LogM(i) = fun_ANMF(R_LogM,x0,s);
        %%%%%% ANMF_LogCC
        Tanmf_LogCC(i) = fun_ANMF(R_LogCC,x0,s);
        %%%%%% ANMF_PCC
        Tanmf_P(i) = fun_ANMF(R_P,x0,s);
        %%%%%% ANMF_PCC
        Tanmf_PCC(i) = fun_ANMF(R_PCC,x0,s);
        %%%%%% ANMF_SFP
        Tanmf_SFP(i) = fun_ANMF(R_SFP,x0,s);  
    end
    TR=sort(Tanmf_R,'descend');
    TCC=sort(Tanmf_CC,'descend');
    TE=sort(Tanmf_E,'descend');
    TECC=sort(Tanmf_ECC,'descend');
    TLogM=sort(Tanmf_LogM,'descend');
    TLogCC=sort(Tanmf_LogCC,'descend');
    TP=sort(Tanmf_P,'descend');
    TPCC=sort(Tanmf_PCC,'descend');
    TSFP = sort(Tanmf_SFP,'descend');

    Th_R(i_Pfa) = (TR(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TR(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_CC(i_Pfa) = (TCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_E(i_Pfa) = (TE(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TE(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_ECC(i_Pfa) = (TECC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TECC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_LogM(i_Pfa) = (TLogM(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TLogM(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_LogCC(i_Pfa) = (TLogCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TLogCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_P(i_Pfa) = (TP(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TP(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_PCC(i_Pfa) = (TPCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TPCC(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    Th_SFP(i_Pfa) = (TSFP(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TSFP(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
end
close(h)
% %%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pd_CC_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_E_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_ECC_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_LogM_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_LogCC_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_P_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_PCC_Mlti_mc = zeros(L_Pfa,length(SNRout));
Pd_SFP_Mlti_mc = zeros(L_Pfa,length(SNRout));
counter_r=0;
counter_cc=0;
counter_e=0;
counter_ecc=0;
counter_logm=0;
counter_logcc=0;
counter_p=0;
counter_pcc=0;
counter_sfp=0;
% alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
alpha=sqrt(SNRnum*mu/(lambda-1));
h = waitbar(0,'Please wait...');
tic
L_SNRout = length(SNRout);
offset = 60000-MonteCarloPd-N+1;
for i_Pfa = 1:L_Pfa %%虚警
   for m=1:L_SNRout %%信噪比
       waitbar(((i_Pfa-1)*L_SNRout+m)/(L_SNRout*L_Pfa),h,sprintf(['检测概率计算: ', num2str(((i_Pfa-1)*L_SNRout+m)/(L_SNRout*L_Pfa)*100),'%%']));
        parfor i=1:MonteCarloPd %%%MC检测概率
    %         waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
            %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%0%
            index_t1 = ceil(rand()*(M-10000))+2000;
            R_KA1 = zeros(N,N);
            R_KA2 = zeros(N,N);
            for ii = 1:before-N+1
                x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
                R_KA1 = R_KA1+fun_NSCMN(x_tt)/(before-N+1);
                R_KA2 = R_KA2+x_tt*x_tt'/(before-N+1);
            end
            Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
            Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
            Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
            x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
            %%%检测信号
            x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
            %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
            R_CC = fun_CC(Train,fun_SCMN(Train),R_KA2);
            R_E = fun_RPowerEMean(Train,1,3);
            R_ECC = fun_PowerCC(Train,R_KA1,1,10);
            R_LogM = fun_RLogEMean(Train,3);
            R_LogCC = fun_LogCC_new(Train,R_KA1,10);
            R_P = fun_RPowerEMean(Train,-1,3);
            R_PCC = fun_PowerCC(Train,R_KA1,-1,9);
            R_SFP = fun_SFP(Train,1);
            
             %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% ANMF_CC
            T_CC = fun_ANMF(R_CC,x0,s);
            %%%%% ANMF_CC
            T_E = fun_ANMF(R_E,x0,s);
            %%%%% ANMF_NSCM
            T_ECC = fun_ANMF(R_ECC,x0,s);
            %%%%%% ANMF_LogCC
            T_LogM = fun_ANMF(R_LogM,x0,s);
            %%%%%% ANMF_LogCC
            T_LogCC = fun_ANMF(R_LogCC,x0,s);
            %%%%%% ANMF_PCC
            T_P = fun_ANMF(R_P,x0,s);
            %%%%%% ANMF_PCC
            T_PCC = fun_ANMF(R_PCC,x0,s);
            %%%%%% ANMF_SFP
            T_SFP = fun_ANMF(R_SFP,x0,s);  
            %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            if T_CC>Th_CC(i_Pfa);               counter_cc=counter_cc+1;    end   
            if T_E>Th_E(i_Pfa);                 counter_e=counter_e+1;    end
            if T_ECC>Th_ECC(i_Pfa);             counter_ecc=counter_ecc+1;    end
            if T_LogM>Th_LogM(i_Pfa);           counter_logm=counter_logm+1;      end
            if T_LogCC>Th_LogCC(i_Pfa);         counter_logcc=counter_logcc+1;      end
            if T_P>Th_P(i_Pfa);                 counter_p=counter_p+1;        end
            if T_PCC>Th_PCC(i_Pfa);             counter_pcc=counter_pcc+1;        end
            if T_SFP>Th_SFP(i_Pfa);             counter_sfp=counter_sfp+1;        end
        end
        Pd_CC_Mlti_mc(i_Pfa,m)=counter_cc/MonteCarloPd;          counter_cc=0;
        Pd_E_Mlti_mc(i_Pfa,m)=counter_e/MonteCarloPd;            counter_e=0;
        Pd_ECC_Mlti_mc(i_Pfa,m)=counter_ecc/MonteCarloPd;        counter_ecc=0;
        Pd_LogM_Mlti_mc(i_Pfa,m)=counter_logm/MonteCarloPd;      counter_logm=0;
        Pd_LogCC_Mlti_mc(i_Pfa,m)=counter_logcc/MonteCarloPd;    counter_logcc=0;
        Pd_P_Mlti_mc(i_Pfa,m)=counter_p/MonteCarloPd;            counter_p=0;
        Pd_PCC_Mlti_mc(i_Pfa,m)=counter_pcc/MonteCarloPd;        counter_pcc=0;
        Pd_SFP_Mlti_mc(i_Pfa,m)=counter_sfp/MonteCarloPd;        counter_sfp=0;
    end
end

close(h)
toc
% figure(2);
% hold on
% %%-5dB
% plot(PFA,Pd_R_Mlti_mc(:,1),'r','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_CC_Mlti_mc(:,1),'g','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_E_Mlti_mc(:,1),'b-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_ECC_Mlti_mc(:,1),'b','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_LogM_Mlti_mc(:,1),'k-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_LogCC_Mlti_mc(:,1),'k','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_P_Mlti_mc(:,1),'c-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_PCC_Mlti_mc(:,1),'c','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_SFP_Mlti_mc(:,1),'m','linewidth',2,'MarkerSize',15)
% h_leg = legend('NMF','CC','E','ECC','LogM','LogCC','P','PCC','SFP');
% xlabel('Pfa','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% %%5dB
% figure(3);
% hold on
% plot(PFA,Pd_R_Mlti_mc(:,2),'r','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_CC_Mlti_mc(:,2),'g','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_E_Mlti_mc(:,2),'b-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_ECC_Mlti_mc(:,2),'b','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_LogM_Mlti_mc(:,2),'k-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_LogCC_Mlti_mc(:,2),'k','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_P_Mlti_mc(:,2),'c-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_PCC_Mlti_mc(:,2),'c','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_SFP_Mlti_mc(:,2),'m','linewidth',2,'MarkerSize',15)
% h_leg = legend('NMF','CC','E','ECC','LogM','LogCC','P','PCC','SFP');
% xlabel('Pfa','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% %%10dB
% figure(4);
% hold on
% plot(PFA,Pd_R_Mlti_mc(:,3),'r','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_CC_Mlti_mc(:,3),'g','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_E_Mlti_mc(:,3),'b-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_ECC_Mlti_mc(:,3),'b','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_LogM_Mlti_mc(:,3),'k-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_LogCC_Mlti_mc(:,3),'k','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_P_Mlti_mc(:,3),'c-.','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_PCC_Mlti_mc(:,3),'c','linewidth',2,'MarkerSize',15)
% plot(PFA,Pd_SFP_Mlti_mc(:,3),'m','linewidth',2,'MarkerSize',15)
% h_leg = legend('NMF','CC','E','ECC','LogM','LogCC','P','PCC','SFP');
% xlabel('Pfa','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
% grid on
% box on
str = ['ROC_IPIX_PFA_',num2str(L),'Second','.mat'];
save (str);
