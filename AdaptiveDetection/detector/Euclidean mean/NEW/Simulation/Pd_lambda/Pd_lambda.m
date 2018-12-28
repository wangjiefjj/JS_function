clc
clear 
close all
%%%%参数设置
n = 1; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 1.5:0.1:3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
rou = 0.95;  %%协方差矩阵生成的迟滞因子
sigma_t =sqrt(0.5);
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
SNRout=10; % 输出SNR
cos2=0.9;
PFA=1e-2;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e3;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
rouR = fun_rho(rou,N,1);
rouR_abs=abs(rouR);
rouR_half=rouR^0.5;
irouR=inv(rouR);
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
Pd_NMF_mc = zeros(1,length(lambda));
Pd_SCM_mc = zeros(1,length(lambda));
Pd_NSCM_mc = zeros(1,length(lambda));
Pd_ECCT_mc = zeros(1,length(lambda));
Pd_ECCS_mc = zeros(1,length(lambda));
Pd_ECCP_mc = zeros(1,length(lambda));
Pd_ML_mc = zeros(1,length(lambda));
Pd_CC_mc = zeros(1,length(lambda));

for i_l = 1:length(lambda)
    i_l
    parfor i = 1:MonteCarloPfa
    %     warning('off')
    %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda(i_l),mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,rouR,lambda(i_l),mu,opt_train); % 接收信号仅包括杂波和噪声
    %     %%%%%%%归一化 
    %     Train = Train/(diag(sqrt(diag(Train'*Train))));
    %     x0 = x0/sqrt(x0'*x0);
        %     %%先验协方差
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
        R_KA =  (rouR).*(t*t');
        R_KA2 =  (tau0*rouR).*(t*t');
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCMN(Train));  
        R_CC = fun_CC(Train,R_SCM,R_KA2);
        R_ML = fun_MLalpha(Train,R_SCM,R_KA2,x0);
        R_ECCT = fun_PowerCC(Train,R_KA,1,4);
        R_ECCS = fun_PowerCC(Train,R_KA,1,8);
        R_ECCP = fun_PowerCC(Train,R_KA,1,7);
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% ANMF_SCM
        Tnmf(i) = fun_ANMF(rouR,x0,s);
        %%%%% ANMF_SCM
        Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,s);
        %%%%% ANMF_ECC
        Tanmf_ECCT(i) = fun_ANMF(R_ECCT,x0,s);
        Tanmf_ECCS(i) = fun_ANMF(R_ECCS,x0,s);
        Tanmf_ECCP(i) = fun_ANMF(R_ECCP,x0,s);
        %%%%%% ANMF_ML
        Tanmf_ML(i) = fun_ANMF(R_ML,x0,s);
        %%%%%% ANMF_CC
        Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
    end
    % close(h)
    toc
    TNMF=sort(Tnmf,'descend');
    TANMF_SCM=sort(Tanmf_SCM,'descend');
    TANMF_ECCT=sort(Tanmf_ECCT,'descend');
    TANMF_ECCS=sort(Tanmf_ECCS,'descend');
    TANMF_ECCP=sort(Tanmf_ECCP,'descend');
    TANMF_ML=sort(Tanmf_ML,'descend');
    TANMF_CC=sort(Tanmf_CC,'descend');

    Th_NMF(i_l) = (TNMF(floor(MonteCarloPfa*PFA-1))+TNMF(floor(MonteCarloPfa*PFA)))/2;
    Th_SCM(i_l) = (TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
    Th_ECCT(i_l)=(TANMF_ECCT(floor(MonteCarloPfa*PFA-1))+TANMF_ECCT(floor(MonteCarloPfa*PFA)))/2;
    Th_ECCS(i_l)=(TANMF_ECCS(floor(MonteCarloPfa*PFA-1))+TANMF_ECCS(floor(MonteCarloPfa*PFA)))/2;
    Th_ECCP(i_l)=(TANMF_ECCP(floor(MonteCarloPfa*PFA-1))+TANMF_ECCP(floor(MonteCarloPfa*PFA)))/2;
    Th_ML(i_l)=(TANMF_ML(floor(MonteCarloPfa*PFA-1))+TANMF_ML(floor(MonteCarloPfa*PFA)))/2;
    Th_CC(i_l) = (TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;

    %%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    counter_nmf=0;
    counter_scm=0;
    counter_nscm=0;
    counter_ecct=0;
    counter_eccs=0;
    counter_eccp=0;
    counter_ml=0;
    counter_cc=0;
    % alpha=sqrt(SNRnum/abs(s'*irouR*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
    if str_train=='g'
        alpha=sqrt(SNRnum);
    elseif str_train=='p'
        alpha=sqrt(SNRnum*mu/(lambda(i_l)-1));     
    end
    parfor i=1:MonteCarloPd
        
      %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,rouR,lambda(i_l),mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,rouR,lambda(i_l),mu,opt_train); % 接收信号仅包括杂波和噪声
        %%先验协方差
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
        R_KA =  (rouR).*(t*t');
        R_KA2 =  (tau0*rouR).*(t*t');
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%
        R_SCM = (fun_SCMN(Train));
        R_ML = fun_MLalpha(Train,R_SCM,R_KA2,x0);    
        R_CC = fun_CC(Train,R_SCM,R_KA2);
        R_ECCT = fun_PowerCC(Train,R_KA,1,4);
        R_ECCS = fun_PowerCC(Train,R_KA,1,8);
        R_ECCP = fun_PowerCC(Train,R_KA,1,7);          
        %%%检测信号
        x0=alpha*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Tnmf = fun_ANMF(rouR,x0,s);
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
        if Tnmf>Th_NMF(i_l);          counter_nmf=counter_nmf+1;        end
        if Tscm>Th_SCM(i_l);          counter_scm=counter_scm+1;        end                
        if Tecct>Th_ECCT(i_l);       counter_ecct=counter_ecct+1;    end
        if Teccs>Th_ECCS(i_l);       counter_eccs=counter_eccs+1;    end
        if Teccp>Th_ECCP(i_l);       counter_eccp=counter_eccp+1;    end
        if Tml>Th_ML(i_l);       counter_ml=counter_ml+1;    end
        if Tcc>Th_CC(i_l);       counter_cc=counter_cc+1;    end
     end
        Pd_NMF_mc(i_l)=counter_nmf/MonteCarloPd;           counter_nmf=0;
        Pd_SCM_mc(i_l)=counter_scm/MonteCarloPd;           counter_scm=0;
        Pd_CC_mc(i_l)=counter_cc/MonteCarloPd;           counter_cc=0;
        Pd_ECCT_mc(i_l)=counter_ecct/MonteCarloPd;        counter_ecct=0;
        Pd_ECCS_mc(i_l)=counter_eccs/MonteCarloPd;        counter_eccs=0;
        Pd_ECCP_mc(i_l)=counter_eccp/MonteCarloPd;        counter_eccp=0;
        Pd_ML_mc(i_l)=counter_ml/MonteCarloPd;        counter_ml=0;
end
figure(1);
hold on
plot(lambda,Pd_NMF_mc,'r','linewidth',2)
plot(lambda,Pd_SCM_mc,'g','linewidth',2)
plot(lambda,Pd_CC_mc,'b','linewidth',2)
plot(lambda,Pd_ML_mc,'c','linewidth',2)
plot(lambda,Pd_ECCT_mc,'k','linewidth',2)
plot(lambda,Pd_ECCS_mc,'K-*','linewidth',2)
plot(lambda,Pd_ECCP_mc,'k-o','linewidth',2)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCT','ANMF with ECCS','ANMF with ECCP');
xlabel('\lambda','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(lambda),max(lambda),0,1])
grid on
% str = [str_train,'_',num2str(n),'N','_s',num2str(sigma_t^2),'.mat'];
% save (str); 

