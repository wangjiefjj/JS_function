%%%%��������
clc
clear 
close all
% Read_Display_Data
Data_process
load(matFile) 
lambda =  2.4072;
mu = 1.3600;
%%%%��������
n = 2; %����������%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=-5:1:25; % ���SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rouR = R_KA;  %%��ʵ���Ӳ�Э����
irouR = inv(rouR);
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��
Zhh = sig;
before = 80+N-1; %%ȥǰ��֡��Ϊ����Э����
tic
% h = waitbar(0,'Please wait...');
parfor i = 1:MonteCarloPfa
    warning('off')
%     warning('query')
%     waitbar(i/MonteCarloPfa,h,sprintf([num2str(i/MonteCarloPfa*100),'%%']));
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    index_t1 = ceil(rand()*(M-10000))+2000;
    %%R_KA
    R_KA = zeros(N,N);
    for ii = 1:before-N+1
        x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
        R_KA = R_KA+x_tt*x_tt'/(before-N+1);
    end
    Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
    Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
    Train = [Train1,Train2];%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = Zhh(index_t1:index_t1+N-1,Range) ; % �����źŽ������Ӳ�������
    %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
%     R_KA = zeros(N,N);
%     for j = 1:before
%         R_KA = R_KA+fun_SCM(Zhh(index_t1-8*j:index_t1-1-8*(j-1), Range))/before;
%     end
%     R_x0 = (fun_SCMN(x0));
     
    R_SCM = (fun_SCMN(Train));
    
    
%     R_NSCM = (fun_NSCMN(Train));

    [R_CC, alpha_cc(i)]= fun_CC(Train,R_SCM,R_KA);
    
    [R_CCIter,alpha_cciter(i)] = fun_CCIter2(Train,R_SCM,R_KA);
    
    [R_ML, alpha_ml(i)] = fun_MLalpha(Train,R_SCM,R_KA,x0);
    
%     R_H = 0.5 * R_KA + 0.5 * R_SCM;
    %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% ANMF_SCM
    Tanmf_SCM(i) = fun_ANMF(R_SCM,x0,s);
    %%%%%% ANMF_NSCM
%     Tanmf_NSCM(i) = fun_ANMF(R_NSCM,x0,s);
    %%%%%% NMF
%     Tnmf(i) = fun_ANMF(rouR,x0,s);
    %%%%%% ANMF_CCIter
    Tanmf_CCIter(i) = fun_ANMF(R_CCIter,x0,s);
    %%%%%% ANMF_ML
    Tanmf_ML(i) = fun_ANMF(R_ML,x0,s);
    %%%%%% ANMF_CC
    Tanmf_CC(i) = fun_ANMF(R_CC,x0,s);
    %%%%%% ANMF_half
%     Tanmf_H(i) = fun_ANMF(R_H,x0,s);
end
toc
TANMF_SCM=sort(Tanmf_SCM,'descend');
% TANMF_NSCM=sort(Tanmf_NSCM,'descend');
% TNMF=sort(Tnmf,'descend');

TANMF_CCIter=sort(Tanmf_CCIter,'descend');
TANMF_ML=sort(Tanmf_ML,'descend');
TANMF_CC=sort(Tanmf_CC,'descend');
% TANMF_H=sort(Tanmf_H,'descend');

Th_SCM = (TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
% Th_NSCM = (TANMF_NSCM(floor(MonteCarloPfa*PFA-1))+TANMF_NSCM(floor(MonteCarloPfa*PFA)))/2;
% Th_NMF = (TNMF(floor(MonteCarloPfa*PFA-1))+TNMF(floor(MonteCarloPfa*PFA)))/2;

Th_CCIter=(TANMF_CCIter(floor(MonteCarloPfa*PFA-1))+TANMF_CCIter(floor(MonteCarloPfa*PFA)))/2;
Th_ML=(TANMF_ML(floor(MonteCarloPfa*PFA-1))+TANMF_ML(floor(MonteCarloPfa*PFA)))/2;
Th_CC = (TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;
% Th_H = (TANMF_H(floor(MonteCarloPfa*PFA-1))+TANMF_H(floor(MonteCarloPfa*PFA)))/2;

m_alpha_ml = mean(alpha_ml);
m_alpha_cc = mean(alpha_cc);
m_alpha_cciter = mean(alpha_cciter);
%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_scm=0;
counter_nscm=0;
counter_nmf=0;
counter_cciter=0;
counter_ml=0;
counter_cc=0;
counter_h=0;

Pd_SCM_mc = zeros(1,length(SNRout));
% Pd_NSCM_mc = zeros(1,length(SNRout));
% Pd_NMF_mc = zeros(1,length(SNRout));
Pd_CCIter_mc = zeros(1,length(SNRout));
Pd_ML_mc = zeros(1,length(SNRout));
Pd_CC_mc = zeros(1,length(SNRout));
% Pd_H_mc = zeros(1,length(SNRout));
% 
% alpha=sqrt(SNRnum/abs(s'*irouR*s)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
alpha=sqrt(SNRnum*mu/(lambda-1));        
h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%0%
        index_t1 = ceil(rand()*(M-10000))+2000;
        %%R_KA
        R_KA = zeros(N,N);
        for ii = 1:before-N+1
            x_tt = Zhh(index_t1-before+ii-1:index_t1-before+ii+N-2,Range);
            R_KA = R_KA+x_tt*x_tt'/(before-N+1);
        end
        Train1 = Zhh(index_t1:index_t1+7,Range-L/2+1:Range-1);
        Train2 = Zhh(index_t1:index_t1+7,Range+1:Range+L/2+1);
        Train = [Train1,Train2];%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = Zhh(index_t1:index_t1+7,Range) ; % �����źŽ������Ӳ�������
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
%         R_KA = zeros(N,N);
%         for j = 1:before
%             R_KA = R_KA+fun_SCM(Zhh(index_t1-8*j:index_t1-1-8*(j-1),Range))/before;
%         end
%         R_x0 = (fun_SCMN(x0));
        
        R_SCM = (fun_SCMN(Train));
    
%         R_NSCM = (fun_NSCMN(Train));
            
        R_CCIter = fun_CCIter2(Train,R_SCM,R_KA);
        
        R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
    
        R_CC = fun_CC(Train,R_SCM,R_KA);
        
%         R_H = 0.5 * R_KA + 0.5 * R_SCM;
        %%%����ź�
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% AMF
        Tscm = fun_ANMF(R_SCM,x0,s);
        %%%%%% ANMF_NSCM
%         Tnscm = fun_ANMF(R_NSCM,x0,s);
        %%%%%% NMF
%         Tnmf = fun_ANMF(rouR,x0,s);
        %%%%%% ANMF_KL
        Tcciter = fun_ANMF(R_CCIter,x0,s);
        %%%%%% ANMF_KL
        Tml = fun_ANMF(R_ML,x0,s);
        %%%%%% ANMF_CC
        Tcc = fun_ANMF(R_CC,x0,s);
        %%%%%% ANMF_CC
%         Thh = fun_ANMF(R_H,x0,s);
        %%%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if Tscm>Th_SCM;          counter_scm=counter_scm+1;        end                
%         if Tnscm>Th_NSCM;       counter_nscm=counter_nscm+1;    end   
%         if Tnmf>Th_NMF;       counter_nmf=counter_nmf+1;    end
        if Tcciter>Th_CCIter;       counter_cciter=counter_cciter+1;    end
        if Tml>Th_ML;       counter_ml=counter_ml+1;    end
        if Tcc>Th_CC;       counter_cc=counter_cc+1;    end
%         if Thh>Th_H;       counter_h=counter_h+1;    end
    end
    Pd_SCM_mc(m)=counter_scm/MonteCarloPd;           counter_scm=0;
%     Pd_NSCM_mc(m)=counter_nscm/MonteCarloPd;        counter_nscm=0;
%     Pd_NMF_mc(m)=counter_nmf/MonteCarloPd;        counter_nmf=0;
    Pd_CC_mc(m)=counter_cc/MonteCarloPd;           counter_cc=0;
    Pd_CCIter_mc(m)=counter_cciter/MonteCarloPd;        counter_cciter=0;
    Pd_ML_mc(m)=counter_ml/MonteCarloPd;        counter_ml=0;
%     Pd_H_mc(m)=counter_h/MonteCarloPd;        counter_h=0;
end
toc
close(h)
figure();
hold on
% plot(SNRout,Pd_NMF_mc,'k','linewidth',1)
plot(SNRout,Pd_SCM_mc,'k-+','linewidth',1)
% plot(SNRout,Pd_NSCM_mc,'k.-','linewidth',1)
plot(SNRout,Pd_CCIter_mc,'k-*','linewidth',1)
plot(SNRout,Pd_ML_mc,'k-o','linewidth',1)
plot(SNRout,Pd_CC_mc,'k-s','linewidth',1)
% plot(SNRout,Pd_H_mc,'k-^','linewidth',1)
h_leg = legend('ANMF with SCM','ANMF with CCIter','ANMF with ML','ANMF with CC');
% h_leg = legend('NMF with correct covariance','NMF with SCM','NMF with NSCM',...
% 'NMF with CCIter','NMF with ML','NMF with CC','NMF with H');

xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
grid on
str = ['R_',num2str(Range),'_','new_CCIter_IPIX_',cdfFile_t,num2str(n),'N','.mat'];
save (str); 

