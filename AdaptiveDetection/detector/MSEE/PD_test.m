clc
clear 
close all
%%参数设置
n = 1; %几倍的样本
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
SCRout=0:1:30; % 输出SCR
SCRnum=10.^(SCRout/10);
CNRout = 30; %%杂噪比
CNRnum=10.^(CNRout/10);
cos2=0.9;
PFA=1e-2;% PFA=1e-4;

MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e3;
L=round(n*N); 
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量

%% 各种协方差
%%杂波协方差1
Rc = CNRnum*fun_rho(rou,N,1,0.1);
R = Rc+eye(N);
%%杂波协方差2
% fc = 0;
% sigmaf = 0.03; %%杂波谱展宽
% rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
% Rc = CNRnum * toeplitz(rc);
% R = Rc+eye(N);
%%%%%%%%%
R_half=R^0.5;
iR=inv(R);
tic
% h = waitbar(1,'Please wait...');
parfor i = 1:MonteCarloPfa

    warning('off')
%%%%%%%%%训练数据产生%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声 
    %%%%%%%%%%%%%%%%%%%%%协方差估计
        R_SCM = (fun_SCMN(Train));
        R_SCM2 = R_SCM + eye(N);
        R_NSCM = (fun_NSCMN(Train));
        R_MESS = fun_MESS(R_SCM);
        R_SFP = fun_SFP(Train,1);
    %%%%%%检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tanmf_R(i) = fun_AMF(R,x0,s);
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_AMF(R_SCM,x0,s);
    %%%%%% ANMF_SCM2
    Tanmf_SCM2(i) = fun_AMF(R_SCM2,x0,s);
    %%%%%% ANMF_NSCM
    Tanmf_NSCM(i) = fun_AMF(R_NSCM,x0,s);
    %%%%%% ANMF_MESS
    Tanmf_MESS(i) = fun_AMF(R_MESS,x0,s);
    %%%%%% ANMF_SFP
    Tanmf_SFP(i) = fun_AMF(R_SFP,x0,s);  
end
% close(h)
toc
TANMF_R=sort(Tanmf_R,'descend');
TANMF_SCM=sort(Tanmf_SCM,'descend');
TANMF_SCM2=sort(Tanmf_SCM2,'descend');
TANMF_NSCM=sort(Tanmf_NSCM,'descend');
TANMF_MESS=sort(Tanmf_MESS,'descend');
TANMF_SFP=sort(Tanmf_SFP,'descend');

Th_R = (TANMF_R(floor(MonteCarloPfa*PFA-1))+TANMF_R(floor(MonteCarloPfa*PFA)))/2;
Th_SCM = (TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
Th_SCM2 = (TANMF_SCM2(floor(MonteCarloPfa*PFA-1))+TANMF_SCM2(floor(MonteCarloPfa*PFA)))/2;
Th_NSCM = (TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;
Th_MESS = (TANMF_MESS(floor(MonteCarloPfa*PFA-1))+TANMF_MESS(floor(MonteCarloPfa*PFA)))/2;
Th_SFP = (TANMF_SFP(floor(MonteCarloPfa*PFA-1))+TANMF_SFP(floor(MonteCarloPfa*PFA)))/2;
%%%%%%%%%%%%%%%%%%%检测概率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_R=0;
counter_SCM=0;
counter_SCM2=0;
counter_NSCM=0;
counter_MESS=0;
counter_SFP=0;

Pd_R_mc = zeros(1,length(SCRout));
Pd_SCM_mc = zeros(1,length(SCRout));
Pd_SCM2_mc = zeros(1,length(SCRout));
Pd_NSCM_mc = zeros(1,length(SCRout));
Pd_MESS_mc = zeros(1,length(SCRout));
Pd_SFP_mc = zeros(1,length(SCRout));
alpha=sqrt(SCRnum/abs(s'*R*s)); % 根据SNR=|alpha|^2*s'*R^(-1)*s求得|alpha|
if str_train=='g'
    alpha=sqrt(SCRnum/abs(s'*Rc*s));
elseif str_train=='p'
    alpha=sqrt(SCRnum*mu/(lambda-1));     
end
% alpha=sqrt(SNRnum*N/abs(s'*s));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SCRout)
    waitbar(m/length(SCRout),h,sprintf([num2str(m/length(SCRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%训练数据产生%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声       
    %%%%%%%%%%%%%%%协方差估计
        R_SCM = (fun_SCMN(Train));
        R_SCM2 = R_SCM + eye(N);
        R_NSCM = (fun_NSCMN(Train));
        R_MESS = fun_MESS(R_SCM);
        R_SFP = fun_SFP(Train,1);
        %检测信号
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
        %检测器%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%OPT
        T_R = fun_AMF(R,x0,s);
        %%%%% ANMF_SCM
        T_SCM = fun_AMF(R_SCM,x0,s);  
        %%%%% ANMF_NSCM
        T_NSCM = fun_AMF(R_NSCM,x0,s); 
        %%%%% ANMF_SCM2
        T_SCM2 = fun_AMF(R_SCM2,x0,s);
        %%%%% ANMF_MESS
        T_MESS = fun_AMF(R_MESS,x0,s);
        %%%%% ANMF_SFP
        T_SFP = fun_AMF(R_SFP,x0,s);  
%         %%判断%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T_R>Th_R;                  counter_R=counter_R+1;          end   
        if T_SCM>Th_SCM;              counter_SCM=counter_SCM+1;      end
        if T_SCM2>Th_SCM2;            counter_SCM2=counter_SCM2+1;    end
        if T_NSCM>Th_NSCM;            counter_NSCM=counter_NSCM+1;    end
        if T_MESS>Th_MESS;            counter_MESS=counter_MESS+1;    end    
        if T_SFP>Th_SFP;              counter_SFP=counter_SFP+1;      end
    end
    Pd_R_mc(m)=counter_R/MonteCarloPd;                counter_R=0;
    Pd_SCM_mc(m)=counter_SCM/MonteCarloPd;            counter_SCM=0;
    Pd_SCM2_mc(m)=counter_SCM2/MonteCarloPd;          counter_SCM2=0;
    Pd_NSCM_mc(m)=counter_NSCM/MonteCarloPd;          counter_NSCM=0;
    Pd_MESS_mc(m)=counter_MESS/MonteCarloPd;          counter_MESS=0;
    Pd_SFP_mc(m)=counter_SFP/MonteCarloPd;            counter_SFP=0;
end
toc
close(h)
figure();
hold on
plot(SCRout,Pd_R_mc,'r.-','linewidth',2)
plot(SCRout,Pd_SCM_mc,'g.-','linewidth',2)
plot(SCRout,Pd_SCM2_mc,'b.-','linewidth',2)
plot(SCRout,Pd_NSCM_mc,'k.-','linewidth',2)
plot(SCRout,Pd_MESS_mc,'c.-','linewidth',2)
plot(SCRout,Pd_SFP_mc,'y.-','linewidth',2)
grid on
box on