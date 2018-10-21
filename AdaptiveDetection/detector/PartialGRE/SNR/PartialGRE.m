clc
clear 
close all
% 
%%
Na=2;     
Np=4;     
N=Na*Np;
optc = 'p';
opt_train = 2;%%1:SIRP,2:部分均匀
L=round(2*N); 
cos2=1;%%%失配
PFA=1e-3;% PFA=1e-4;
%%各种比
SNRout = 0:30; % 输出信噪比SNR
CNRout = 15; %杂噪比
JNRout = 15; %干噪比
SNRnum=10.^(SNRout/10);
CNRnum=10.^(CNRout/10);
JNRnum=10.^(JNRout/10);
%%虚警率和MC次数
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
%%信号
fs=[0.3];%,0.3175
theta = ones(length(fs),1);
nn=(0:N-1)';
H = exp(-1i*2*pi*nn*fs); %%%%%% 导向矢量
vt = H*theta;

%% 各种协方差
%%杂波协方差1
fc = 0;
sigmaf = 0.03; %%杂波谱展宽
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc = CNRnum * toeplitz(rc);
%%杂波协方差2
% fc = 0;
% rho = 0.9;
% Rc = CNRnum * fun_rho(rho,N,1,fc);
%%干扰协方差
fj = -0.2;
jam = exp(-1i*2*pi*(nn-N/2)*fj);
Rj = JNRnum * jam*jam';
%%检测单元协方差test
Rt = eye(N) + Rc + Rj;
%%参考单元协方差secondary 
Rs = (eye(N) + 10*Rc);
%% 失配导向矢量设置
iRs = inv(Rs);
[UU,SS,VV]=svd(iRs*vt);
s_v=UU(:,2); %%%%%% 与vt在白化空间正交，即：vt^H*iR*s_v==0
weight=linspace(0,1,300);
for i=1:length(weight)
    s_tmpt=weight(i)*vt+(1-weight(i))*s_v;
    cos2_tmpt(i)=abs(s_tmpt'*iRs*vt).^2/abs(s_tmpt'*iRs*s_tmpt*vt'*iRs*vt);
end
[Min, Index]=min(abs(cos2-cos2_tmpt));
Weight=weight(Index);
s_real=Weight*vt+(1-Weight)*s_v;
%% 门限
Tparglrt = zeros(1,MonteCarloPfa);
Tsglrt = zeros(1,MonteCarloPfa);
Tanmf = zeros(1,MonteCarloPfa);
tic    
% h = waitbar(0,'Please wait...'); 
parfor i=1:MonteCarloPfa
    warning off
    %% 数据产生
    Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
    x = fun_TrainData(optc,N,1,Rt ,3,1,opt_train);
    %%
    R_SCM = fun_SCMN(Train);
    R_NSCM = fun_NSCMN(Train);
    %% 检测器
    Tparglrt(i)=fun_Partial_GER_GLRT(Train,x,H);                   %%%%%% KGLRT
    Tsglrt(i)=fun_SGLRT(Train,x,H);                   %%%%%% KGLRT
    Tanmf(i)=fun_ANMF(R_NSCM,x,H)
end
toc

TPARGLRT=sort(Tparglrt,'descend');
TSGLRT=sort(Tsglrt,'descend');
TANMF=sort(Tanmf,'descend');


Th_PARGLRT = (TPARGLRT(floor(MonteCarloPfa*PFA-1))+TPARGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_SGLRT = (TSGLRT(floor(MonteCarloPfa*PFA-1))+TSGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_ANMF = (TANMF(floor(MonteCarloPfa*PFA-1))+TANMF(floor(MonteCarloPfa*PFA)))/2;
toc
a=0;b=2*pi;
%% 开始检测
counter_sglrt=0;
counter_parglrt=0;
counter_anmf=0;
%% 信噪比
alpha=sqrt(SNRnum*N/abs(vt'*vt));
%% 检测概率 

Pd_SGLRT_mc = zeros(1,length(SNRout));
Pd_PARGLRT_mc = zeros(1,length(SNRout));
Pd_ANMF_mc = zeros(1,length(SNRout));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd
        warning off
        %% 数据产生
        Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
        x = fun_TrainData(optc,N,1,Rt,3,1,opt_train);
        x=alpha(m)*s_real+x;%+pp;    %%%%%%%  閲嶈  %%%%%%%%%%%%%
       %%
        R_SCM = fun_SCMN(Train);
        R_NSCM = fun_NSCMN(Train);
       %% 检测器
        Tparglrt=fun_Partial_GER_GLRT(Train,x,H);                   %%%%%% SGLRT
        Tsglrt = fun_SGLRT(Train,x,H);               %%%SGLRT
        Tanmf=fun_ANMF(R_NSCM,x,H)
       %% 比较
        if Tparglrt>Th_PARGLRT;         counter_parglrt=counter_parglrt+1;          end 
        if Tsglrt>Th_SGLRT;         counter_sglrt=counter_sglrt+1;          end 
        if Tanmf>Th_ANMF;         counter_anmf=counter_anmf+1;          end 
    end
    Pd_PARGLRT_mc(m)=counter_parglrt/MonteCarloPd;          counter_parglrt=0;
    Pd_SGLRT_mc(m)=counter_sglrt/MonteCarloPd;          counter_sglrt=0;
    Pd_ANMF_mc(m)=counter_anmf/MonteCarloPd;          counter_anmf=0;
end
close(h)
toc

figure(2);
hold on
plot(SNRout,Pd_PARGLRT_mc,'r-.x','linewidth',2,'MarkerSize',10)
plot(SNRout,Pd_SGLRT_mc,'k-.o','linewidth',2,'MarkerSize',2)
plot(SNRout,Pd_ANMF_mc,'b-.x','linewidth',2,'MarkerSize',10)
legend('PartialGLRT','SGLRT','ANMF')
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on
box on