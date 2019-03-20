%%%GLC-GLRT的ROC曲线
clc
clear 
close all
% %%%%参数设置
n = 0.5; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_tt = 0.9;
sigma_t =sqrt(sigma_tt);
%%%Pd_CLGLRT_2Kmu1lambda3s0.1o1_p：2K：训练单元数目，mu，lambda，s：失配向量方差，
%%o1:opt=1，p：IG纹理复合高斯
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=[-5,5,10]; % 输出SNR
cos2=0.9;
PFA=[1e-4,1e-3:1e-2:1e-1+1e-2];% PFA=1e-4;[1e-4,1e-3:1e-2:1e-1+1e-2];
SNRnum=10.^(SNRout/10);
MonteCarloPfa=round(1./PFA*100);
L_Pfa = length(MonteCarloPfa);
MonteCarloPd=1e4;
rou = 0.90;  %%协方差矩阵生成的迟滞因子
rouR = zeros(N,N);  %%真实的杂波协方差
L=round(n*N); 
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% 系统导向矢量
for i=1:N
    for j=1:N
        rouR(i,j)=rou^abs(i-j);%*exp(1j*2*pi*abs(i-j)*theta_sig);
    end
end
irouR=inv(rouR);
rouR_abs=abs(rouR);
R_KA = zeros(size(rouR));
tic
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
for i  = 1:10000
    R_KA = R_KA+rouR.*(t*t')/10000;
end
iR_KA = inv(R_KA);
toc
% R_KA_inv = inv(R_KA);
rouR_half=rouR^0.5;
%%%%导向矢量设置
% [UU,SS,VV]=svd(irouR*s);
% s_v=UU(:,2); %%%%%% 与vt在白化空间正交，即：s^H*iR*s_v==0
% weight=linspace(0,1,300);
% for i=1:length(weight)
%     s_tmpt=weight(i)*s+(1-weight(i))*s_v;
%     cos2_tmpt(i)=abs(s_tmpt'*irouR*s).^2/abs(s_tmpt'*irouR*s_tmpt*s'*irouR*s);
% end
% [Min, Index]=min(abs(cos2-cos2_tmpt));
% Weight=weight(Index);
% s_real=Weight*s+(1-Weight)*s_v;
% figure;plot(abs(s_real))
% figure; plot(weight,cos2_tmpt);
% %%%%%正式开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%门限计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    parfor i = 1:MonteCarloPfa(i_Pfa)
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量    
        R_KACC = (R).*(t*t');
        R_SCM = fun_SCMN(Train);
        R_CC = fun_CC(Train,R_SCM,R_KACC);
        R_E = fun_RPowerEMean(Train,1,10);
        R_ECC = fun_PowerCC(Train,R_KA,1,10);
        R_LogM = fun_RLogEMean(Train,10);
        R_LogCC = fun_LogCC_new(Train,R_KA,10);
        R_P = fun_RPowerEMean(Train,-1,10);
        R_PCC = fun_PowerCC(Train,R_KA,-1,10);
        R_SFP = fun_SFP(Train,1);
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%OPT
        Tanmf_R(i) = fun_ANMF(rouR,x0,s);
        %%%%%% ANMF_CC
        Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
        %%%%%% ANMF_CC
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
Pd_R_Mlti_mc = zeros(L_Pfa,length(SNRout));
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
% alpha=sqrt(SNRnum/abs(s_real'*irouR*s_real)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
alpha=sqrt(SNRnum*mu/(lambda-1));
h = waitbar(0,'Please wait...');
tic
L_SNRout = length(SNRout);
for i_Pfa = 1:L_Pfa %%虚警
   for m=1:L_SNRout %%信噪比
       waitbar(((i_Pfa-1)*L_SNRout+m)/(L_SNRout*L_Pfa),h,sprintf(['检测概率计算: ', num2str(((i_Pfa-1)*L_SNRout+m)/(L_SNRout*L_Pfa)*100),'%%']));
        parfor i=1:MonteCarloPd %%%MC检测概率
    %         waitbar(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd,h,sprintf([num2str(((m-1)*MonteCarloPd+i)/length(SNRout)/MonteCarloPd*100),'%%']));
            %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
            Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
            x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
           %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
            t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量    
            R_KACC = (R).*(t*t');
            R_SCM = fun_SCMN(Train);
            R_CC = fun_CC(Train,R_SCM,R_KA);
            R_E = fun_RPowerEMean(Train,1,10);
            R_ECC = fun_PowerCC(Train,R_KA,1,10);
            R_LogM = fun_RLogEMean(Train,10);
            R_LogCC = fun_LogCC_new(Train,R_KA,10);
            R_P = fun_RPowerEMean(Train,-1,10);
            R_PCC = fun_PowerCC(Train,R_KA,-1,10);
            R_SFP = fun_SFP(Train,1);
            x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
             %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%OPT
            T_R = fun_ANMF(rouR,x0,s);
            %%%%%% ANMF_CC
            T_CC = fun_ANMF(R_CC,x0,s);
            %%%%%% ANMF_CC
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
            if T_R>Th_R(i_Pfa);                 counter_r=counter_r+1;        end                  
            if T_CC>Th_CC(i_Pfa);               counter_cc=counter_cc+1;    end   
            if T_E>Th_E(i_Pfa);                 counter_e=counter_e+1;    end
            if T_ECC>Th_ECC(i_Pfa);             counter_ecc=counter_ecc+1;    end
            if T_LogM>Th_LogM(i_Pfa);           counter_logm=counter_logm+1;      end
            if T_LogCC>Th_LogCC(i_Pfa);         counter_logcc=counter_logcc+1;      end
            if T_P>Th_P(i_Pfa);                 counter_p=counter_p+1;        end
            if T_PCC>Th_PCC(i_Pfa);             counter_pcc=counter_pcc+1;        end
            if T_SFP>Th_SFP(i_Pfa);             counter_sfp=counter_sfp+1;        end
        end
        Pd_R_Mlti_mc(i_Pfa,m)=counter_r/MonteCarloPd;            counter_r=0;
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
figure(2);
hold on
%%-5dB
plot(PFA,Pd_R_Mlti_mc(:,1),'r','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_CC_Mlti_mc(:,1),'g','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_E_Mlti_mc(:,1),'b-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_ECC_Mlti_mc(:,1),'b','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_LogM_Mlti_mc(:,1),'k-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_LogCC_Mlti_mc(:,1),'k','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_P_Mlti_mc(:,1),'c-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_PCC_Mlti_mc(:,1),'c','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_SFP_Mlti_mc(:,1),'m','linewidth',2,'MarkerSize',15)
h_leg = legend('NMF','CC','E','ECC','LogM','LogCC','P','PCC','SFP');
xlabel('Pfa','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
box on
%%5dB
figure(3);
hold on
plot(PFA,Pd_R_Mlti_mc(:,2),'r','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_CC_Mlti_mc(:,2),'g','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_E_Mlti_mc(:,2),'b-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_ECC_Mlti_mc(:,2),'b','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_LogM_Mlti_mc(:,2),'k-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_LogCC_Mlti_mc(:,2),'k','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_P_Mlti_mc(:,2),'c-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_PCC_Mlti_mc(:,2),'c','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_SFP_Mlti_mc(:,2),'m','linewidth',2,'MarkerSize',15)
h_leg = legend('NMF','CC','E','ECC','LogM','LogCC','P','PCC','SFP');
xlabel('Pfa','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
box on
%%10dB
figure(4);
hold on
plot(PFA,Pd_R_Mlti_mc(:,3),'r','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_CC_Mlti_mc(:,3),'g','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_E_Mlti_mc(:,3),'b-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_ECC_Mlti_mc(:,3),'b','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_LogM_Mlti_mc(:,3),'k-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_LogCC_Mlti_mc(:,3),'k','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_P_Mlti_mc(:,3),'c-.','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_PCC_Mlti_mc(:,3),'c','linewidth',2,'MarkerSize',15)
plot(PFA,Pd_SFP_Mlti_mc(:,3),'m','linewidth',2,'MarkerSize',15)
h_leg = legend('NMF','CC','E','ECC','LogM','LogCC','P','PCC','SFP');
xlabel('Pfa','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
box on
str = ['ROC_',num2str(L),'Second','_s',num2str(sigma_t),'_',str_train,'.mat'];
save (str);
