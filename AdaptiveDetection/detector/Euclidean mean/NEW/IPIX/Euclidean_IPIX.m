%%%%参数设置
clc
clear 
close all
% Read_Display_Data
Data_process
load(matFile) 
lambda =  2.4072; %17,1.1967, 19980204_224024 %8,2.4072, 19980223_170435
mu = 1.3600;      %17,1.3180, 19980204_224024 %8,1.3600, 19980223_170435
%%%%参数设置
n = 2;      %几倍的样本%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=-5:1:25; % 输出SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rouR = R_KA;  %%真实的杂波协方差
irouR = inv(rouR);
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
Zhh = sig;
before = 100+N-1; %%去前几帧作为先验协方差
tic
% h = waitbar(0,'Please wait...');
parfor i = 1:MonteCarloPfa
    warning('off')
%     warning('query')
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
    index_t1 = ceil(rand()*(M-10000))+2000;
    %%R_KA
    R_KA = zeros(N,N);
    for ii = 1:before-N+1
        x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
        x_tt = x_tt/norm(x_tt);
        R_KA = R_KA+x_tt*x_tt'/(before-N+1);
    end
    Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
    Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
    Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % 接收信号仅包括杂波和噪声
    %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCMN(Train));  
    R_CC = fun_CC(Train,R_SCM,R_KA);
    R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
    R_ECCT = fun_PowerCC(Train,R_KA,1,4);
    R_ECCS = fun_PowerCC(Train,R_KA,1,8);
    R_ECCP = fun_PowerCC(Train,R_KA,1,7);
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,s);
    %%%%%% ANMF_CCIter
    Tanmf_ECCT(i) = fun_ANMF(R_ECCT,x0,s);
    Tanmf_ECCS(i) = fun_ANMF(R_ECCS,x0,s);
    Tanmf_ECCP(i) = fun_ANMF(R_ECCP,x0,s);
    %%%%%% ANMF_ML
    Tanmf_ML(i) = fun_ANMF(R_ML,x0,s);
    %%%%%% ANMF_CC
    Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
end
toc
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_ECCT=sort(Tanmf_ECCT,'descend');
TANMF_ECCS=sort(Tanmf_ECCS,'descend');
TANMF_ECCP=sort(Tanmf_ECCP,'descend');
TANMF_ML=sort(Tanmf_ML,'descend');
TANMF_CC=sort(Tanmf_CC,'descend');

Th_SCM = (TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_ECCT=(TANMF_ECCT(floor(MonteCarloPfa*PFA-1))+TANMF_ECCT(floor(MonteCarloPfa*PFA)))/2;
Th_ECCS=(TANMF_ECCS(floor(MonteCarloPfa*PFA-1))+TANMF_ECCS(floor(MonteCarloPfa*PFA)))/2;
Th_ECCP=(TANMF_ECCP(floor(MonteCarloPfa*PFA-1))+TANMF_ECCP(floor(MonteCarloPfa*PFA)))/2;
Th_ML=(TANMF_ML(floor(MonteCarloPfa*PFA-1))+TANMF_ML(floor(MonteCarloPfa*PFA)))/2;
Th_CC = (TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;

%%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_scm=0;
counter_nscm=0;
counter_ecct=0;
counter_eccs=0;
counter_eccp=0;
counter_ml=0;
counter_cc=0;

Pd_SCM_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_ECCT_mc = zeros(1,length(SNRout));
Pd_ECCS_mc = zeros(1,length(SNRout));
Pd_ECCP_mc = zeros(1,length(SNRout));
Pd_ML_mc = zeros(1,length(SNRout));
Pd_CC_mc = zeros(1,length(SNRout));
Pd_H_mc = zeros(1,length(SNRout));
% 
% alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
alpha=sqrt(SNRnum*mu/(lambda-1));        
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%0%
        index_t1 = ceil(rand()*(M-10000))+2000;
        %%R_KA
        R_KA = zeros(N,N);
        for ii = 1:before-N+1
            x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
            R_KA = R_KA+x_tt*x_tt'/(before-N+1);
        end
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
        Train = [Train1,Train2];%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = Zhh(index_t1:index_t1+7,Range) ; % 接收信号仅包括杂波和噪声
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
%%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
        R_SCM = (fun_SCMN(Train));    
        R_ECCT = fun_PowerCC(Train,R_KA,1,4);
        R_ECCS = fun_PowerCC(Train,R_KA,1,8);
        R_ECCP = fun_PowerCC(Train,R_KA,1,7);
        R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);    
        R_CC = fun_CC(Train,R_SCM,R_KA);  
        %%%检测信号
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% AMF
        Tscm = fun_ANMF(R_SCM,x0,s);
        %%%%%% ANMF_CC
        Tecct = fun_ANMF(R_ECCT,x0,s);
        Teccs = fun_ANMF(R_ECCS,x0,s);
        Teccp = fun_ANMF(R_ECCP,x0,s);
        %%%%%% ANMF_ML
        Tml = fun_ANMF(R_ML,x0,s);
        %%%%%% ANMF_CC
        Tcc = fun_ANMF(R_CC,x0,s);
        %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if Tscm>Th_SCM;          counter_scm=counter_scm+1;        end                
        if Tecct>Th_ECCT;       counter_ecct=counter_ecct+1;    end
        if Teccs>Th_ECCS;       counter_eccs=counter_eccs+1;    end
        if Teccp>Th_ECCP;       counter_eccp=counter_eccp+1;    end
        if Tml>Th_ML;       counter_ml=counter_ml+1;    end
        if Tcc>Th_CC;       counter_cc=counter_cc+1;    end
    end
    Pd_SCM_mc(m)=counter_scm/MonteCarloPd;           counter_scm=0;
    Pd_CC_mc(m)=counter_cc/MonteCarloPd;           counter_cc=0;
    Pd_ECCT_mc(m)=counter_ecct/MonteCarloPd;        counter_ecct=0;
    Pd_ECCS_mc(m)=counter_eccs/MonteCarloPd;        counter_eccs=0;
    Pd_ECCP_mc(m)=counter_eccp/MonteCarloPd;        counter_eccp=0;
    Pd_ML_mc(m)=counter_ml/MonteCarloPd;        counter_ml=0;
end
toc
close(h)
figure(2);
hold on
plot(SNRout,Pd_SCM_mc,'g','linewidth',2)
plot(SNRout,Pd_CC_mc,'b','linewidth',2)
plot(SNRout,Pd_ML_mc,'c','linewidth',2)
plot(SNRout,Pd_ECCT_mc,'k','linewidth',2)
plot(SNRout,Pd_ECCS_mc,'K-*','linewidth',2)
plot(SNRout,Pd_ECCP_mc,'k-o','linewidth',2)
h_leg = legend('ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with T-KA','ANMF with S-KA','ANMF with P-KA');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
str = ['IPIX_',cdfFile_t,num2str(n),'N','.mat'];
save (str); 


