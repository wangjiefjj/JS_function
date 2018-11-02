%%%检测器在不同检测环境下的性能
clc
clear 
close all
%%%%参数设置
n = 2; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
%%%%假设参数设置
SNRout=0:1:35; % 输出SNR
SNRnum=10.^(SNRout/10);
CNRout=30;
CNRnum=10.^(CNRout/10);
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
theta_sig = 0.3;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量

%%仿真参数设置
PFA=1e-3;% PFA=1e-4;
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;  
L=round(n*N); 
rou = 0.95;  %%协方差矩阵生成的迟滞因子
%%杂波协方差
fc=0;
sigmaf = 0.03; %%杂波谱展宽
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc1 = CNRnum*toeplitz(rc);
Rc1 = Rc1 + eye(N);%+ eye(N)
Rc2 = Rc1;
% rouR = fun_rho(rou,N,1,0.1);%%真实的杂波协方差
% irouR = inv(rouR);
tic
% h = waitbar(1,'Please wait...');
parfor i = 1:MonteCarloPfa
%     warning('off')
%%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,Rc1);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    [x0,tau0] = fun_TrainData(str_train,N,1,Rc2); % 接收信号仅包括杂波和噪声
    %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
    R_SCM = (fun_SCM(Train));  
    R_NSCM = (fun_NSCMN(Train));  
    R_AML = fun_AML(Train);
    %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% KGLRT
    Tamf(i) = fun_AMF(R_SCM,x0,s);
    %%%%% ANMF_NSCM
    Tanmf_NSCM(i) = fun_ANMF(R_SCM,x0,s);
    %%%%% ANMF_AML
    Tanmf_AML(i) = fun_ANMF(R_AML,x0,s);
end
% close(h)
toc
TAMF=sort(Tamf,'descend');
TANMF_NSCM=sort(Tanmf_NSCM,'descend');
TANMF_AML=sort(Tanmf_AML,'descend');


Th_AMF = (TAMF(floor(MonteCarloPfa*PFA-1))+TAMF(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM=(TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;
Th_AML=(TANMF_AML(floor(MonteCarloPfa*PFA-1))+TANMF_AML(floor(MonteCarloPfa*PFA)))/2;

%%%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_amf=0;
counter_nscm=0;
counter_aml=0;

Pd_AMF_mc = zeros(1,length(SNRout));
Pd_NSCM_mc = zeros(1,length(SNRout));
Pd_AML_mc = zeros(1,length(SNRout));
alpha=sqrt(SNRnum); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,Rc1);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,Rc2); % 接收信号仅包括杂波和噪声
        %%%%协方差估计%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCM(Train));  
        R_NSCM = (fun_NSCMN(Train));  
        R_AML = fun_AML(Train);      
        %%%检测信号
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% KGLRT
        Th_amf = fun_AMF(R_SCM,x0,s);
        %%%%% ANMF_NSCM
        Th_nscm = fun_ANMF(R_SCM,x0,s);
        %%%%% ANMF_AML
        Thaml = fun_ANMF(R_AML,x0,s);
        %%%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        if Th_amf>Th_AMF;          counter_amf=counter_amf+1;        end                
        if Th_nscm>Th_NSCM;       counter_nscm=counter_nscm+1;    end
        if Thaml>Th_AML;       counter_aml=counter_aml+1;    end
    end   
    Pd_AMF_mc(m)=counter_amf/MonteCarloPd;           counter_amf=0;
    Pd_NSCM_mc(m)=counter_nscm/MonteCarloPd;        counter_nscm=0;
    Pd_AML_mc(m)=counter_aml/MonteCarloPd;        counter_aml=0;
end
toc
close(h)
figure(1);
hold on
plot(SNRout,Pd_AMF_mc,'g','linewidth',2)
plot(SNRout,Pd_NSCM_mc,'K-*','linewidth',2)
plot(SNRout,Pd_AML_mc,'k','linewidth',2)
h_leg = legend('KGLRT','NSCM','AML');
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
% str = ['g','_',num2str(n),'N','.mat'];
% save (str); 

