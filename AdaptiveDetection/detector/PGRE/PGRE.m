clc
clear 
close all
% 
%%
Na=2;     
Np=4;     
N=Na*Np;
optc = 'g';
opt_train = 1;%%1:SIRP,2:部分均匀
L=round(1*N); 
cos2=1;%%%失配
PFA=1e-3;% PFA=1e-4;
%%各种比
SNRout = 0:1:30; % 输出信噪比SNR
CNRout = 15; %杂噪比
JNRout = 15; %干噪比
SNRnum=10.^(SNRout/10);
CNRnum=10.^(CNRout/10);
JNRnum=10.^(JNRout/10);
%%虚警率和MC次数
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
%%信号
fs=[0.3,0.3175];
theta = ones(length(fs),1);
nn=(0:N-1)';
H = exp(-1i*2*pi*nn*fs); %%%%%% 导向矢量
vt = H*theta;

% [UU,SS,VV]=svd(iR*H);
% vt_v=UU(:,end); %导向矢量和杂波的正交分量vt^H*iR*vt_v==0,
% q = 2*vt_v;%%GRE中的q
% R0 = R + q*q';
%% 各种协方差
%%杂波协方差
fc = 0;
sigmaf = 0.03; %%杂波谱展宽
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc = CNRnum * toeplitz(rc);
%%干扰协方差
fj = -0.2;
jam = exp(-1i*2*pi*nn*fs);
Rj = JNRnum * jam*jam';
%%检测单元协方差test
Rt = eye(N) + Rc + Rj;
%%参考单元协方差secondary 
Rs = eye(N) + Rc;
%% 门限
Tpglrt = zeros(1,MonteCarloPfa);
Tprao = zeros(1,MonteCarloPfa);
Tpwald = zeros(1,MonteCarloPfa);
tic    
% h = waitbar(0,'Please wait...'); 
for i=1:MonteCarloPfa
    warning off
    %% 数据产生
    Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
    x = fun_TrainData(optc,N,1,Rt ,3,1,opt_train);
    %% 检测器
    Tpglrt(i)=fun_P_GRE_GLRT(Train,x,H);                   %%%%%% KGLRT
    Tprao(i)=fun_P_GRE_Rao(Train,x,H);                     %%%%%% RAO
    Tpwald(i)=fun_P_GRE_Wald(Train,x,H);                   %%%%%% Wald
 
end
toc
TPGLRT=sort(Tpglrt,'descend');
TPRAO=sort(Tprao,'descend');
TPWALD=sort(Tpwald,'descend');

Th_PGLRT = (TPGLRT(floor(MonteCarloPfa*PFA-1))+TPGLRT(floor(MonteCarloPfa*PFA)))/2;
Th_PRAO = (TPRAO(floor(MonteCarloPfa*PFA-1))+TPRAO(floor(MonteCarloPfa*PFA)))/2;
Th_PWALD=(TPWALD(floor(MonteCarloPfa*PFA-1))+TPWALD(floor(MonteCarloPfa*PFA)))/2;
toc
a=0;b=2*pi;
%% 开始检测
tic
counter_pglrt=0;
counter_prao=0;
counter_pwald=0;

%% 信噪比
alpha=sqrt(SNRnum*N/abs(vt'*vt));
%% 检测概率 
Pd_PGLRT_mc = zeros(1,length(SNRout));
Pd_PRAO_mc = zeros(1,length(SNRout));
Pd_PWALD_mc = zeros(1,length(SNRout));
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
    parfor i=1:MonteCarloPd
        warning off
        %% 数据产生
        Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
        x = fun_TrainData(optc,N,1,Rt,3,1,opt_train);
        x=alpha(m)*vt+x;%+pp;    %%%%%%%  閲嶈  %%%%%%%%%%%%%
        %% 检测器
        Tpglrt=fun_P_GRE_GLRT(Train,x,H);                   %%%%%% KGLRT
        Tprao=fun_P_GRE_Rao(Train,x,H);                     %%%%%% RAO
        Tpwald=fun_P_GRE_Wald(Train,x,H);                   %%%%%% Wald
        %%
        if Tpglrt>Th_PGLRT;         counter_pglrt=counter_pglrt+1;          end 
        if Tprao>Th_PRAO;         counter_prao=counter_prao+1;          end 
        if Tpwald>Th_PWALD;         counter_pwald=counter_pwald+1;          end 
    end
    Pd_PGLRT_mc(m)=counter_pglrt/MonteCarloPd;          counter_pglrt=0;
    Pd_PRAO_mc(m)=counter_prao/MonteCarloPd;          counter_prao=0;
    Pd_PWALD_mc(m)=counter_pwald/MonteCarloPd;          counter_pwald=0;
end
close(h)
toc

figure(2);
hold on
plot(SNRout,Pd_PGLRT_mc,'k-o','linewidth',2,'MarkerSize',2)
plot(SNRout,Pd_PRAO_mc,'r-x','linewidth',2,'MarkerSize',10)
plot(SNRout,Pd_PWALD_mc,'b-x','linewidth',2,'MarkerSize',10)
legend('PGLRT','PRAO','PWALD')
xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on