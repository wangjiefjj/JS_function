clc
clear 
close all
% ����ȶԼ�����ܵ�Ӱ��
%%
Na=2;     
Np=4;     
N=Na*Np;
optc = 'g';
opt_train = 1;%%1:SIRP,2:���־���
L=round(2*N); 
cos2=1;%%%ʧ��
PFA=1e-3;% PFA=1e-4;
%%���ֱ�
SNRout = 25; % ��������SNR
CNRout = 15; %�����
JNRout = 0:40; %�����
SNRnum=10.^(SNRout/10);
CNRnum=10.^(CNRout/10);
JNRnum=10.^(JNRout/10);
%%�龯�ʺ�MC����
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
%%�ź�
fs=[0.3,0.3175];
theta = ones(length(fs),1);
nn=(0:N-1)';
H = exp(-1i*2*pi*nn*fs)/sqrt(N); %%%%%% ����ʸ��
vt = H*theta;
%% �����
alpha=sqrt(SNRnum*N/abs(vt'*vt));
%% ����Э����
%%�Ӳ�Э����
fc = 0;
sigmaf = 0.03; %%�Ӳ���չ��
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
Rc = CNRnum * toeplitz(rc);
%%�ο���ԪЭ����secondary 
Rs = eye(N) + Rc;
%% ������ 
Pd_PGLRT_mc = zeros(1,length(JNRnum));
Pd_PRAO_mc = zeros(1,length(JNRnum));
Pd_PWALD_mc = zeros(1,length(JNRnum));
Th_PGLRT = zeros(1,length(JNRnum));
Th_PRAO = zeros(1,length(JNRnum));
Th_PWALD = zeros(1,length(JNRnum));

Pd_SGLRT_mc = zeros(1,length(JNRnum));
Pd_SRAO_mc = zeros(1,length(JNRnum));
Pd_SWALD_mc = zeros(1,length(JNRnum));
Th_SGLRT = zeros(1,length(JNRnum));
Th_SRAO = zeros(1,length(JNRnum));
Th_SWALD = zeros(1,length(JNRnum));
%% ����ȶԼ�����ܵ�Ӱ�����
h = waitbar(0,'Please wait...'); 
for i_JNR = 1:length(JNRnum)
    waitbar(i_JNR/length(JNRnum),h,sprintf([num2str(i_JNR/length(JNRnum)*100),'%%']));
    %%����Э����
    fj = -0.2;
    jam = exp(-1i*2*pi*nn*fj);
    Rj = JNRnum(i_JNR) * (jam*jam');
    %%��ⵥԪЭ����test
    Rt = eye(N) + Rc + Rj;
    %% ����
    Tpglrt = zeros(1,MonteCarloPfa);
    Tprao = zeros(1,MonteCarloPfa);
    Tpwald = zeros(1,MonteCarloPfa);  
    Tsglrt = zeros(1,MonteCarloPfa);
    Tsrao = zeros(1,MonteCarloPfa);
    Tswald = zeros(1,MonteCarloPfa);
    parfor i=1:MonteCarloPfa
        warning off
        %% ���ݲ���
        Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
        x = fun_TrainData(optc,N,1,Rt ,3,1,opt_train);
        %% �����
        Tpglrt(i)=fun_P_GRE_GLRT(Train,x,H);                   %%%%%% KGLRT
        Tprao(i)=fun_P_GRE_Rao(Train,x,H);                     %%%%%% RAO
        Tpwald(i)=fun_P_GRE_Wald(Train,x,H);                   %%%%%% Wald
        Tsglrt(i)=fun_SGLRT(Train,x,H);                       %%%%%% KGLRT
        Tsrao(i)=fun_SRAO(Train,x,H);                         %%%%%% RAO
        Tswald(i)=fun_SWALD(Train,x,H);                       %%%%%% Wald
    end
    TPGLRT=sort(Tpglrt,'descend');
    TPRAO=sort(Tprao,'descend');
    TPWALD=sort(Tpwald,'descend');
    
    TSGLRT=sort(Tsglrt,'descend');
    TSRAO=sort(Tsrao,'descend');
    TSWALD=sort(Tswald,'descend');
    
    Th_PGLRT(i_JNR) = (TPGLRT(floor(MonteCarloPfa*PFA-1))+TPGLRT(floor(MonteCarloPfa*PFA)))/2;
    Th_PRAO(i_JNR) = (TPRAO(floor(MonteCarloPfa*PFA-1))+TPRAO(floor(MonteCarloPfa*PFA)))/2;
    Th_PWALD(i_JNR)=(TPWALD(floor(MonteCarloPfa*PFA-1))+TPWALD(floor(MonteCarloPfa*PFA)))/2;
    
    Th_SGLRT(i_JNR) = (TSGLRT(floor(MonteCarloPfa*PFA-1))+TSGLRT(floor(MonteCarloPfa*PFA)))/2;
    Th_SRAO(i_JNR) = (TSRAO(floor(MonteCarloPfa*PFA-1))+TSRAO(floor(MonteCarloPfa*PFA)))/2;
    Th_SWALD(i_JNR)=(TSWALD(floor(MonteCarloPfa*PFA-1))+TSWALD(floor(MonteCarloPfa*PFA)))/2;
    %% ��ʼ���
    counter_pglrt=0;
    counter_prao=0;
    counter_pwald=0;
    
    counter_sglrt=0;
    counter_srao=0;
    counter_swald=0;
    parfor i=1:MonteCarloPd
    warning off
        %% ���ݲ���
        Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
        x = fun_TrainData(optc,N,1,Rt,3,1,opt_train);
        x=alpha*vt+x;%+pp;    %%%%%%%  重要  %%%%%%%%%%%%%
       %%
        Tpglrt=fun_P_GRE_GLRT(Train,x,H);                   %%%%%% PGLRT
        Tprao = fun_P_GRE_Rao(Train,x,H);               %%%PRAO
        Tpwald = fun_P_GRE_Wald(Train,x,H);               %%%PWald
        
        Tsglrt=fun_SGLRT(Train,x,H);                   %%%%%% SGLRT
        Tsrao = fun_SRAO(Train,x,H);               %%%SRAO
        Tswald = fun_SWALD(Train,x,H);               %%%SWald
       %%
        if Tpglrt>Th_PGLRT(i_JNR);         counter_pglrt=counter_pglrt+1;          end 
        if Tprao>Th_PRAO(i_JNR);           counter_prao=counter_prao+1;            end 
        if Tpwald>Th_PWALD(i_JNR);         counter_pwald=counter_pwald+1;          end 
        
        if Tsglrt>Th_SGLRT(i_JNR);         counter_sglrt=counter_sglrt+1;          end 
        if Tsrao>Th_SRAO(i_JNR);           counter_srao=counter_srao+1;            end 
        if Tswald>Th_SWALD(i_JNR);         counter_swald=counter_swald+1;          end 
        
    end
    Pd_PGLRT_mc(i_JNR)=counter_pglrt/MonteCarloPd;          counter_pglrt=0;
    Pd_PRAO_mc(i_JNR)=counter_prao/MonteCarloPd;            counter_prao=0;
    Pd_PWALD_mc(i_JNR)=counter_pwald/MonteCarloPd;          counter_pwald=0;
    
    Pd_SGLRT_mc(i_JNR)=counter_sglrt/MonteCarloPd;          counter_sglrt=0;
    Pd_SRAO_mc(i_JNR)=counter_srao/MonteCarloPd;            counter_srao=0;
    Pd_SWALD_mc(i_JNR)=counter_swald/MonteCarloPd;          counter_swald=0;
end
close(h)

figure(2);
hold on
plot(JNRout,Pd_PGLRT_mc,'k-o','linewidth',2,'MarkerSize',2)
plot(JNRout,Pd_PRAO_mc,'r-x','linewidth',2,'MarkerSize',10)
plot(JNRout,Pd_PWALD_mc,'b-x','linewidth',2,'MarkerSize',10)

plot(JNRout,Pd_SGLRT_mc,'k-.o','linewidth',2,'MarkerSize',2)
plot(JNRout,Pd_SRAO_mc,'r-.x','linewidth',2,'MarkerSize',10)
plot(JNRout,Pd_SWALD_mc,'b-.x','linewidth',2,'MarkerSize',10)
legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD')
xlabel('JNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
grid on
box on