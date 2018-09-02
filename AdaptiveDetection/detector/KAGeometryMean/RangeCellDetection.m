%%每个距离单元的检测情况,RD图
clc
clear 
close all
%%%%参数设置
% n = 1.25; %几倍的样本

str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
rou = 0.90;  %%协方差矩阵生成的迟滞因子
sigma_t = 0.9;
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
Sample_num = 10;
Sample_num_half = Sample_num/2;
N = Na*Np;
SNRout= 10; % 输出SNR
cos2=0.9;
PFA=1e-2;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e3;
rouR = zeros(N,N);  %%真实的杂波协方差
L=40;%round(n*N);%%距离单元个数
LL = 3*L;
nn = 0:N-1;
rouR = fun_rho(rou,N,1,0);
rouR_abs=abs(rouR);
rouR_half=rouR^0.5;
irouR=inv(rouR);
t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
R_KA = zeros(size(rouR));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%失配向量
    R_KA = R_KA + rouR.*(t*t')/1000;
end
if str_train=='g'
    alpha=sqrt(SNRnum);
elseif str_train=='p'
    alpha=sqrt(SNRnum*mu/(lambda-1));     
end
%%原始杂波加目标图
MC = 10
theta_sig = 0.2%-0.5:0.1:0.5;
theta_t = 0.2;
index1 = find( abs(theta_sig - theta_t)<1e-5);
TSV = zeros(length(theta_sig),L);
TSV_SCM = zeros(length(theta_sig),L);
TSV_CC = zeros(length(theta_sig),L);
TSV_ML = zeros(length(theta_sig),L);
TSV_NSCM = zeros(length(theta_sig),L);
TSV_ECC = zeros(length(theta_sig),L);
TSV_LogM = zeros(length(theta_sig),L);
TSV_LogCC = zeros(length(theta_sig),L);
TSV_P = zeros(length(theta_sig),L);
TSV_PCC = zeros(length(theta_sig),L);
h = waitbar(0,'Please wait...');
for i_s = 1:length(theta_sig)
    for mc = 1:MC
        waitbar(((i_s-1)*MC+mc)/(length(theta_sig)*MC),h,sprintf([num2str(((i_s-1)*MC+mc)/(length(theta_sig)*MC)*100),'%%']));
        Data = fun_TrainData(str_train,N,LL,rouR,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        s = exp(-1i*2*pi*nn*theta_t).'/sqrt(N); %%%%%% 目标导向矢量
        p = exp(-1i*2*pi*nn*theta_sig(i_s)).'/sqrt(N); %%%%%% 系统导向矢量
        x0 = Data(:,L+L/2-1); % 接收信号仅包括杂波和噪声
        x0= alpha*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        Data(:,L+L/2-1) = x0;
        for i = L:2*L-1
            Train = [Data(:,i-Sample_num_half:i-1),Data(:,i+1:i+Sample_num_half)];
            CUT = Data(:,i);
            %%协方差估计
            R_SCM = fun_SCMN(Train);
            R_CC = fun_CC(Train,R_SCM,R_KA);
            R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
            R_NSCM = fun_NSCMN(Train);
            R_ECC = fun_CCIter2(Train,R_NSCM,R_KA);
            R_LogM = fun_RLogEMean(Train,4);
            R_LogCC = fun_LogCC_new(Train,R_KA,4);
            R_P = fun_RPowerEMean(Train,-1,4);
            R_PCC = fun_PowerCC(Train,R_KA,-1,4);  
            %%检验统计量
            %%%%
            TSV(i_s,i-L+1) = TSV(i_s,i-L+1)+fun_ANMF(rouR,CUT,p)/MC;
            %%%%
            t_SCM = fun_ANMF(R_SCM,CUT,p);
            if t_SCM>1
                t_SCM = 0;
            end           
            TSV_SCM(i_s,i-L+1) = TSV_SCM(i_s,i-L+1)+t_SCM/MC;
            %%%%
            TSV_CC(i_s,i-L+1) = TSV_CC(i_s,i-L+1)+fun_ANMF(R_CC,CUT,p)/MC;
            %%%%
            t_ML = fun_ANMF(R_ML,CUT,p);
            if t_ML>1
                t_ML = 0;
            end
            TSV_ML(i_s,i-L+1) = TSV_ML(i_s,i-L+1)+t_ML/MC;
            %%%%
            t_NSCM = fun_ANMF(R_NSCM,CUT,p);
            if t_NSCM>1
                t_NSCM = 0;
            end
            TSV_NSCM(i_s,i-L+1) = TSV_NSCM(i_s,i-L+1)+t_NSCM/MC;%TestStatisticValue
            %%%%
            TSV_ECC(i_s,i-L+1) = TSV_ECC(i_s,i-L+1)+fun_ANMF(R_ECC,CUT,p)/MC;
            TSV_LogM(i_s,i-L+1) = TSV_LogM(i_s,i-L+1)+fun_ANMF(R_LogM,CUT,p)/MC;
            TSV_LogCC(i_s,i-L+1) = TSV_LogCC(i_s,i-L+1)+fun_ANMF(R_LogCC,CUT,p)/MC;
            TSV_P(i_s,i-L+1) = TSV_P(i_s,i-L+1)+fun_ANMF(R_P,CUT,p)/MC;
            TSV_PCC(i_s,i-L+1) = TSV_PCC(i_s,i-L+1)+fun_ANMF(R_PCC,CUT,p)/MC; 
        end
    end
end
close(h)
[X,Y] = meshgrid(1:L,theta_sig);
% if length(theta_sig) ==1
    figure
    hold on
    plot(TSV_SCM(index1,:),'b-.')%/max(TSV_SCM(index1,:))
    plot(TSV_CC(index1,:),'b-.*')%/max(TSV_CC(index1,:))
    plot(TSV_ML(index1,:),'b-.>')%/max(TSV_ML(index1,:))
    plot(TSV_NSCM(index1,:),'k-.')%/max(TSV_NSCM(index1,:))
    plot(TSV_ECC(index1,:),'k-.*')%/max(TSV_ECC(index1,:))
    plot(TSV_LogM(index1,:),'r.-')%/max(TSV_LogM(index1,:))
    plot(TSV_LogCC(index1,:),'r-.*')%/max(TSV_LogCC(index1,:))
    plot(TSV_P(index1,:),'g.-')%/max(TSV_P(index1,:))
    plot(TSV_PCC(index1,:),'g-.*')%/max(TSV_PCC(index1,:))
    plot(TSV(index1,:),'m-*')%/max(TSV(index1,:))
    grid on
    h_leg = legend('SCM','CC','ML','NSCM','ECC','LogM','LogCC','P','PCC','OPT');
% else
%     figure
%     mesh(X,Y,TSV_NSCM/max(max(TSV_NSCM)))
%     xlabel('Range Cell')
%     ylabel('Normalized Doppler')
%     zlabel('Normalized Statistics')
%     figure
%     mesh(X,Y,TSV_LogM/max(max(TSV_LogM)))
%     xlabel('Range Cell')
%     ylabel('Normalized Doppler')
%     zlabel('Normalized Statistics')
%     figure
%     mesh(X,Y,TSV_LogCC/max(max(TSV_LogCC)))
%     xlabel('Range Cell')
%     ylabel('Normalized Doppler')
%     zlabel('Normalized Statistics')
%     figure
%     mesh(X,Y,TSV_ECC/max(max(TSV_ECC)))
%     xlabel('Range Cell')
%     ylabel('Normalized Doppler')
%     zlabel('Normalized Statistics')
%     figure
%     mesh(X,Y,TSV_CC/max(max(TSV_CC))) 
%     xlabel('Range Cell')
%     ylabel('Normalized Doppler')
%     zlabel('Normalized Statistics')
% end
str = ['TSV_',num2str(Sample_num),'Second',num2str(SNRout),'SNR','_s',num2str(sigma_t),'_',str_train,'.mat'];
save (str); 