clc
clear 
close all
%%参数设置
n = 0.5; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
rou = 0.90;  %%协方差矩阵生成的迟滞因子
sigma_t =0.1;
%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=10; % 输出SNR
CNR = 30; %%杂噪比
cos2=0.9;
SNRnum=10.^(SNRout/10);
PFA=1e-3;% PFA=1e-4;
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e3;
L=round(n*N); 
nn = 0:N-1;
fd=0.2;
s = exp(-1i*2*pi*nn*fd).'/sqrt(N); %%%%%% 系统导向矢量
fc = -0.5:0.05:0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%门限计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_fc = 1:length(fc)
    i_fc/length(fc)
    rouR = fun_rho(rou,N,1,fc(i_fc)); %%真实的杂波协方差
    rouR_half=rouR^0.5;
    irouR=inv(rouR);
    R_KA1 = zeros(size(rouR));
    for i = 1:1000
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
        R_KA1 = R_KA1 + rouR.*(t*t')/1000;
    end
    Tanmf_R = zeros(1,MonteCarloPfa);
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
    %%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    %     Train = awgn(Train,CNR);
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
    %     x0 = awgn(Train,CNR);
    % % % %     RKA
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
        R_KA2 = rouR.*(t*t');    
    % % % %     协方差估计
        R_CC = fun_CC(Train,fun_SCMN(Train),R_KA2);
        R_E = fun_RPowerEMean(Train,1,3);
        R_ECC = fun_PowerCC(Train,R_KA1,1,4);
        R_LogM = fun_RLogEMean(Train,3);
        R_LogCC = fun_LogCC_new(Train,R_KA1,4);
        R_P = fun_RPowerEMean(Train,-1,3);
        R_PCC = fun_PowerCC(Train,R_KA1,-1,4);
        R_SFP = fun_SFP(Train,1);
    %     %%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Tanmf_R(i) = fun_ANMF(rouR,x0,s);
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
    TANMF_R=sort(Tanmf_R,'descend');
    TANMF_CC=sort(Tanmf_CC,'descend');
    TANMF_E=sort(Tanmf_E,'descend');
    TANMF_ECC=sort(Tanmf_ECC,'descend');
    TANMF_LogM=sort(Tanmf_LogM,'descend');
    TANMF_LogCC=sort(Tanmf_LogCC,'descend');
    TANMF_P=sort(Tanmf_P,'descend');
    TANMF_PCC=sort(Tanmf_PCC,'descend');
    TANMF_SFP=sort(Tanmf_SFP,'descend');

    Th_R(i_fc) = (TANMF_R(floor(MonteCarloPfa*PFA-1))+TANMF_R(floor(MonteCarloPfa*PFA)))/2;
    Th_CC(i_fc) = (TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;
    Th_E(i_fc) = (TANMF_E(floor(MonteCarloPfa*PFA-1))+TANMF_E(floor(MonteCarloPfa*PFA)))/2;
    Th_ECC(i_fc) = (TANMF_ECC(floor(MonteCarloPfa*PFA-1))+TANMF_ECC(floor(MonteCarloPfa*PFA)))/2;
    Th_LogM(i_fc) = (TANMF_LogM(floor(MonteCarloPfa*PFA-1))+TANMF_LogM(floor(MonteCarloPfa*PFA)))/2;
    Th_LogCC(i_fc) = (TANMF_LogCC(floor(MonteCarloPfa*PFA-1))+TANMF_LogCC(floor(MonteCarloPfa*PFA)))/2;
    Th_P(i_fc) = (TANMF_P(floor(MonteCarloPfa*PFA-1))+TANMF_P(floor(MonteCarloPfa*PFA)))/2;
    Th_PCC(i_fc) = (TANMF_PCC(floor(MonteCarloPfa*PFA-1))+TANMF_PCC(floor(MonteCarloPfa*PFA)))/2;
    Th_SFP(i_fc) = (TANMF_SFP(floor(MonteCarloPfa*PFA-1))+TANMF_SFP(floor(MonteCarloPfa*PFA)))/2;
end

% load Th_4Second_s0.1_PFA3_p.mat
%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0,'Please wait...');
%     alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
if str_train=='g'
   alpha=sqrt(SNRnum/abs(s'*irouR*s));
elseif str_train=='p'
    alpha=sqrt(SNRnum*mu/(lambda-1));     
end
Pd_R_mc = zeros(1,length(fc));
Pd_CC_mc = zeros(1,length(fc));
Pd_E_mc = zeros(1,length(fc));
Pd_ECC_mc = zeros(1,length(fc));
Pd_LogM_mc = zeros(1,length(fc));
Pd_LogCC_mc = zeros(1,length(fc));
Pd_P_mc = zeros(1,length(fc));
Pd_PCC_mc = zeros(1,length(fc));
Pd_SFP_mc = zeros(1,length(fc));
tic
for i_t = 1:length(fc)
    waitbar(i_t/length(fc),h,sprintf([num2str(i_t/length(fc)*100),'%%']));
    rouR = fun_rho(rou,N,1,fc(i_t)); %%真实的杂波协方差
    rouR_half=rouR^0.5;
    irouR=inv(rouR);
    R_KA1 = zeros(size(rouR));
    for i = 1:1000
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
        R_KA1 = R_KA1 + rouR.*(t*t')/1000;
    end
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
%         Train = awgn(Train,CNR);
        x0 = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
%         x0 = awgn(x0,CNR);
% % %       RKA
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
        R_KA2 = rouR.*(t*t');        
% %         协方差估计
        R_CC = fun_CC(Train,fun_SCMN(Train),R_KA2);
        R_E = fun_RPowerEMean(Train,1,3);
        R_ECC = fun_PowerCC(Train,R_KA1,1,4);
        R_LogM = fun_RLogEMean(Train,3);
        R_LogCC = fun_LogCC_new(Train,R_KA1,4);
        R_P = fun_RPowerEMean(Train,-1,3);
        R_PCC = fun_PowerCC(Train,R_KA1,-1,4);
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
%         %%% ANMF_LogM
        T_LogM = fun_ANMF(R_LogM,x0,s);
        if T_LogM>1
            T_LogM = 1;
        end
        %%%% ANMF_LogCC
        T_LogCC = fun_ANMF(R_LogCC,x0,s);
        %%%% ANMF_P
        T_P = fun_ANMF(R_P,x0,s);
        %%%%% ANMF_PCC
        T_PCC = fun_ANMF(R_PCC,x0,s);
        %%%%% ANMF_SFP
        T_SFP = fun_ANMF(R_SFP,x0,s);  
%         %%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T_R>Th_R(i_t);                counter_R=counter_R+1;    end   
        if T_CC>Th_CC(i_t);              counter_CC=counter_CC+1;    end
        if T_E>Th_E(i_t);            counter_E=counter_E+1;    end
        if T_ECC>Th_ECC(i_t);            counter_ECC=counter_ECC+1;    end
        if T_LogM>Th_LogM(i_t);          counter_LogM=counter_LogM+1;    end
        if T_LogCC>Th_LogCC(i_t);        counter_LogCC=counter_LogCC+1;    end
        if T_P>Th_P(i_t);                counter_P=counter_P+1;    end
        if T_PCC>Th_PCC(i_t);            counter_PCC=counter_PCC+1;    end
        if T_SFP>Th_SFP(i_t);            counter_SFP=counter_SFP+1;    end
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
plot(fc,Pd_R_mc,'c.-','linewidth',1)
plot(fc,Pd_CC_mc,'b-*','linewidth',1)
plot(fc,Pd_E_mc,'k-.*','linewidth',1)
plot(fc,Pd_ECC_mc,'k-*','linewidth',1)
plot(fc,Pd_LogM_mc,'r.-','linewidth',1)
plot(fc,Pd_LogCC_mc,'r-*','linewidth',1)
plot(fc,Pd_P_mc,'g.-','linewidth',1)
plot(fc,Pd_PCC_mc,'g-*','linewidth',1)
plot(fc,Pd_SFP_mc,'c-*','linewidth',1)

h_leg = legend('NMF','ANMF with CC',...
    'ANMF with E','ANMF with ECC','ANMF with LogM','ANMF with LogCC',...
    'ANMF with P','ANMF with PCC','SFP');

% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
% grid on
% grid minor
% str = ['PdFd_2_',num2str(L),'Second','_',str_train,'.mat'];
str = ['PdFc_',num2str(L),'Second','_s',num2str(sigma_t),'_',str_train,'.mat'];
% % % % % str = ['Pd_E',num2str(L),'Second','_',str_train,'.mat'];
save (str); 
