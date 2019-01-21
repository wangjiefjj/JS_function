clc
clear 
close all
% 
%%
Na=2;     
Np=8;     
N=Na*Np;
optc = 'p';
opt_train = 1;%%1:SIRP,2:部分均匀
lambda = 2;
mu = 1;
L=round(2*N); 
cos2=1;%%%失配
PFA=[1e-3:1e-3:1e-2,1e-1-1e-2:1e-2:1e-1];% PFA=1e-4;
%%各种比
SNRout = 0:1:30; % 输出信噪比SNR
CNRout = 25; %杂噪比
JNRout = 15; %干噪比
SNRnum=10.^(SNRout/10);
CNRnum=10.^(CNRout/10);
JNRnum=10.^(JNRout/10);
%%虚警率和MC次数
MonteCarloPfa=round(1./PFA*100);
MonteCarloPFA=1e4;
%%信号
fd=0.2;
theta = ones(length(fd),1);
nn=(0:N-1)';
vt = exp(1i*2*pi*nn*fd)/sqrt(N); %%%%%% 导向矢量

%% 各种协方差
%%杂波协方差
% fc = 0.15;
% sigmaf = 0.03; %%杂波谱展宽
% rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
% Rc = CNRnum * toeplitz(rc);
Rc = CNRnum * fun_rho(0.9,N,1,0.1);  
%%干扰协方差
fj = -0.2;
jam = exp(-1i*2*pi*nn*fj);
Rj = JNRnum * jam*jam';
%%检测单元协方差test
Rt = eye(N) + Rc;
%%参考单元协方差secondary 
Rs = eye(N) + Rc;
%% 门限
PFA_ANMF_mc = zeros(1,length(PFA));
PFA_MTI_mc = zeros(1,length(PFA));
h = waitbar(0,'Please wait...'); 
tic
for i_Pfa = 1:length(PFA)
    waitbar(i_Pfa/length(PFA),h,sprintf([num2str(i_Pfa/length(PFA)*100),'%%']));
    Tamf = zeros(1,MonteCarloPfa(i_Pfa));
    parfor i=1:MonteCarloPfa(i_Pfa)
        warning off
        %% 数据产生
        Train = fun_TrainData(optc,N,L,Rs,lambda,mu,opt_train);
        x = fun_TrainData(optc,N,1,Rt,lambda,mu,opt_train);
        RSCM = fun_NSCMN(Train);
        %% 检测器
        Tamf(i)=fun_ANMF(RSCM,x,vt);                 %%%%%% AMF
    end
    TAMF=sort(Tamf,'descend');
    Th_ANMF(i_Pfa) = (TAMF(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa)-1))+TAMF(floor(MonteCarloPfa(i_Pfa)*PFA(i_Pfa))))/2;
    %% 实际虚警概率 
    counter_amf=0;
    counter_mti=0;
    %% CFAR,18个单元的CA-CFAR
    Reference = min(18,L-6);
    guard = 2;%%保护单元
    L_t = Reference/2+guard;
    cfar = ones(Reference+2*guard+1,1)/Reference; 
    CUT = Reference/2+3;
    cfar(CUT-2:CUT+2)=0;
    alpha=(Reference)*(PFA(i_Pfa)^(-1/(Reference))-1);%求alpha
%%
%% 实际虚警率
parfor i=1:MonteCarloPfa(i_Pfa)
    Train = fun_TrainData(optc,N,L,Rs,lambda,mu,opt_train);
    x = fun_TrainData(optc,N,1,Rt,lambda,mu,opt_train);
    %% 自适应检测
    RSCM = fun_NSCMN(Train);
    %%检测器
    Tamf=fun_ANMF(RSCM,x,vt);                   %%%%%% ANMF
    %% MTI
    Pc = [Train(:,1:L/2),x,Train(:,L/2+1:end)];%%包含目标的脉压结果
    MTI = fun_MTI(Pc);
    [Lfft,~]=size(MTI);
    target_d = round(Lfft/2+fd/0.5*Lfft*0.5)+1;%%目标所在多普勒通道
    %%CA-CFAR 
    Tmti = MTI(target_d,L/2+1);
    Mean = sum(MTI(target_d,L/2-L_t:L/2+L_t).'.*cfar);
    Th_MTI=alpha*Mean;%%求门限
    %
    if Tamf>Th_ANMF(i_Pfa);         counter_amf=counter_amf+1;          end 
    if Tmti>Th_MTI;         counter_mti=counter_mti+1;          end 
end
PFA_ANMF_mc(i_Pfa)=counter_amf/MonteCarloPfa(i_Pfa);          counter_amf=0;
PFA_MTI_mc(i_Pfa)=counter_mti/MonteCarloPfa(i_Pfa);          counter_mti=0;
end
toc


%%
figure(2);
hold on
plot(1:length(PFA),PFA,'g','linewidth',2,'MarkerSize',2)
plot(1:length(PFA),PFA_ANMF_mc,'k-o','linewidth',2,'MarkerSize',2)
plot(1:length(PFA),PFA_MTI_mc,'r-x','linewidth',2,'MarkerSize',10)
legend('设定值','ANMF实际虚警率','MTI-FFT-CFAR实际虚警率')
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on
box on