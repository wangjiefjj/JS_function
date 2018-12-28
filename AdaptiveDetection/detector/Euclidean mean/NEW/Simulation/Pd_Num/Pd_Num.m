clc
clear 
close all
%%%%��������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�������ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
tau_m = mu/(lambda-1);
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG��������ͬ
rou = 0.95;  %%Э����������ɵĳ�������
sigma_t =sqrt(0.5);
%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=10; % ���SNR
cos2=0.9;
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
rouR = zeros(N,N);  %%��ʵ���Ӳ�Э����
L=4:16; 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��
rouR = fun_rho(rou,N,1);
rouR_abs=abs(rouR);
rouR_half=rouR^0.5;
irouR=inv(rouR);
t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
tic
% h = waitbar(1,'Please wait...');
Pd_NMF_mc = zeros(1,length(L));
Pd_SCM_mc = zeros(1,length(L));
Pd_NSCM_mc = zeros(1,length(L));
Pd_ECCT_mc = zeros(1,length(L));
Pd_ECCS_mc = zeros(1,length(L));
Pd_ECCP_mc = zeros(1,length(L));
Pd_ML_mc = zeros(1,length(L));
Pd_CC_mc = zeros(1,length(L));
h = waitbar(0,'Please wait...');
for i_L=1:length(L)
    waitbar(i_L/length(L),h,sprintf([num2str(i_L/length(L)*100),'%%']));
    
    parfor i = 1:MonteCarloPfa
%     warning('off')
%%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L(i_L),rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        [x0,tau0] = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
    %     %%%%%%%��һ�� 
    %     Train = Train/(diag(sqrt(diag(Train'*Train))));
    %     x0 = x0/sqrt(x0'*x0);
        %     %%����Э����
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
        R_KA =  (rouR).*(t*t');
        R_KA2 =  (tau0*rouR).*(t*t');
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCMN(Train));  
        R_CC = fun_CC(Train,R_SCM,R_KA2);
        R_ML = fun_MLalpha(Train,R_SCM,R_KA2,x0);
        R_ECCT = fun_PowerCC(Train,R_KA,1,4);
        R_ECCS = fun_PowerCC(Train,R_KA,1,8);
        R_ECCP = fun_PowerCC(Train,R_KA,1,7);
        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    TNMF=sort(Tnmf,'descend');
    TANMF_SCM=sort(Tanmf_SCM,'descend');
    TANMF_ECCT=sort(Tanmf_ECCT,'descend');
    TANMF_ECCS=sort(Tanmf_ECCS,'descend');
    TANMF_ECCP=sort(Tanmf_ECCP,'descend');
    TANMF_ML=sort(Tanmf_ML,'descend');
    TANMF_CC=sort(Tanmf_CC,'descend');

    Th_NMF(i_L) = (TNMF(floor(MonteCarloPfa*PFA-1))+TNMF(floor(MonteCarloPfa*PFA)))/2;
    Th_SCM(i_L) = (TANMF_SCM(floor(MonteCarloPfa*PFA-1))+TANMF_SCM(floor(MonteCarloPfa*PFA)))/2;
    Th_ECCT(i_L)=(TANMF_ECCT(floor(MonteCarloPfa*PFA-1))+TANMF_ECCT(floor(MonteCarloPfa*PFA)))/2;
    Th_ECCS(i_L)=(TANMF_ECCS(floor(MonteCarloPfa*PFA-1))+TANMF_ECCS(floor(MonteCarloPfa*PFA)))/2;
    Th_ECCP(i_L)=(TANMF_ECCP(floor(MonteCarloPfa*PFA-1))+TANMF_ECCP(floor(MonteCarloPfa*PFA)))/2;
    Th_ML(i_L) =(TANMF_ML(floor(MonteCarloPfa*PFA-1))+TANMF_ML(floor(MonteCarloPfa*PFA)))/2;
    Th_CC(i_L) = (TANMF_CC(floor(MonteCarloPfa*PFA-1))+TANMF_CC(floor(MonteCarloPfa*PFA)))/2;
    %%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    counter_nmf=0;
    counter_scm=0;
    counter_nscm=0;
    counter_ecct=0;
    counter_eccs=0;
    counter_eccp=0;
    counter_ml=0;
    counter_cc=0;
    if str_train=='g'
        alpha=sqrt(SNRnum);
    elseif str_train=='p'
        alpha=sqrt(SNRnum*mu/(lambda-1));     
    end
    parfor i=1:MonteCarloPd
        %%%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L(i_L),rouR,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        [x0,tau0] = fun_TrainData(str_train,N,1,rouR,lambda,mu,opt_train); % �����źŽ������Ӳ�������
        %%����Э����
        t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
        R_KA =  (rouR).*(t*t');
        R_KA2 =  (tau0*rouR).*(t*t');
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%
        R_SCM = (fun_SCMN(Train));    
        R_ECCT = fun_PowerCC(Train,R_KA,1,4);
        R_ECCS = fun_PowerCC(Train,R_KA,1,8);
        R_ECCP = fun_PowerCC(Train,R_KA,1,7);
        R_ML = fun_MLalpha(Train,R_SCM,R_KA2,x0);    
        R_CC = fun_CC(Train,R_SCM,R_KA2);        
        %%%����ź�
        x0=alpha*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %%%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if Tnmf>Th_NMF(i_L);          counter_nmf=counter_nmf+1;        end
        if Tscm>Th_SCM(i_L);          counter_scm=counter_scm+1;        end                
        if Tecct>Th_ECCT(i_L);       counter_ecct=counter_ecct+1;    end
        if Teccs>Th_ECCS(i_L);       counter_eccs=counter_eccs+1;    end
        if Teccp>Th_ECCP(i_L);       counter_eccp=counter_eccp+1;    end
        if Tml>Th_ML(i_L);       counter_ml=counter_ml+1;    end
        if Tcc>Th_CC(i_L);       counter_cc=counter_cc+1;    end
    end
    Pd_NMF_mc(i_L)=counter_nmf/MonteCarloPd;           counter_nmf=0;
    Pd_SCM_mc(i_L)=counter_scm/MonteCarloPd;           counter_scm=0;
    Pd_CC_mc(i_L)=counter_cc/MonteCarloPd;           counter_cc=0;
    Pd_ECCT_mc(i_L)=counter_ecct/MonteCarloPd;        counter_ecct=0;
    Pd_ECCS_mc(i_L)=counter_eccs/MonteCarloPd;        counter_eccs=0;
    Pd_ECCP_mc(i_L)=counter_eccp/MonteCarloPd;        counter_eccp=0;
    Pd_ML_mc(i_L)=counter_ml/MonteCarloPd;        counter_ml=0;
end

figure(1);
hold on
plot(L,Pd_NMF_mc,'r','linewidth',2)
plot(L,Pd_SCM_mc,'g','linewidth',2)
plot(L,Pd_CC_mc,'b','linewidth',2)
plot(L,Pd_ML_mc,'c','linewidth',2)
plot(L,Pd_ECCT_mc,'k','linewidth',2)
plot(L,Pd_ECCS_mc,'K-*','linewidth',2)
plot(L,Pd_ECCP_mc,'k-o','linewidth',2)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCT','ANMF with ECCS','ANMF with ECCP');
xlabel('Number of secondary data','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
set(h_leg,'Location','SouthEast')
axis([min(L),max(L),0,1])
grid on
str = [str_train,'_Num_','_s',num2str(sigma_t^2),'.mat'];
save (str); 
