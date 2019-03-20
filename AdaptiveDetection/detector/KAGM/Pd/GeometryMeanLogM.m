clc
clear 
close all
%%��������
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
tau_m = mu/(lambda-1);
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
rou = 0.90;  %%Э����������ɵĳ�������
sigma_t =0.9;
%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
SNRout=-5:1:25; % ���SNR
CNRout = 30; %%�����
CNRnum=10.^(CNRout/10);
cos2=0.9;
PFA=1e-1;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e2;
L=round(n*N); 
theta_sig = 0.2;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% ϵͳ����ʸ��
%%�Ӳ�Э����
fc=0;
Rc = fun_rho(rou,N,1,fc);
R = Rc;
% R = CNRnum*Rc+eye(N);
R_KA = zeros(size(R));
for i = 1:10000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
    R_KA = R_KA + R.*(t*t')/10000;
end
tic
% h = waitbar(1,'Please wait...');
parfor i = 1:MonteCarloPfa

    warning('off')
%%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
%     Train = awgn(Train,CNR);
    [x0,tau0] = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); % �����źŽ������Ӳ�������  
% % % %     Э�������
    R_LogM = fun_RLogEMean(Train,10);
    R_LogCC = fun_LogCC_new(Train,R_KA,10);
%     %%�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ANMF_LogM
    Tanmf_LogM(i) = fun_ANMF(R_LogM,x0,s);
    if Tanmf_LogM(i)>1
        Tanmf_LogM(i) = 0;
    end
    %%%% ANMF_LogCC
    Tanmf_LogCC(i) = fun_ANMF(R_LogCC,x0,s);
end
% close(h)
toc
TANMF_LogM=sort(Tanmf_LogM,'descend');
TANMF_LogCC=sort(Tanmf_LogCC,'descend');

Th_LogM = (TANMF_LogM(floor(MonteCarloPfa*PFA-1))+TANMF_LogM(floor(MonteCarloPfa*PFA)))/2;
Th_LogCC = (TANMF_LogCC(floor(MonteCarloPfa*PFA-1))+TANMF_LogCC(floor(MonteCarloPfa*PFA)))/2;
%%
%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_LogM=0;
counter_LogCC=0;

Pd_LogM_mc = zeros(1,length(SNRout));
Pd_LogCC_mc = zeros(1,length(SNRout));
% alpha=sqrt(SNRnum/abs(s'*irouR*s)); % ����SNR=|alpha|^2*s'*R^(-1)*s���|alpha|
if str_train=='g'
    alpha=sqrt(SNRnum/abs(s'*irouR*s));
elseif str_train=='p'
    alpha=sqrt(SNRnum*mu/(lambda-1));     
end

h = waitbar(0,'Please wait...');
tic
for m=1:length(SNRout)
    waitbar(m/length(SNRout),h,sprintf([num2str(m/length(SNRout)*100),'%%']));
%     MonteCarloPd = 1e4;
    parfor i=1:MonteCarloPd
        %%%%%%%%%ѵ�����ݲ���%%%%%%%%%%%%%%
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
%         Train = awgn(Train,CNR);
        [x0,tau0] = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); % �����źŽ������Ӳ�������
%         x0 = awgn(x0,CNR);
%         Э�������
        R_LogM = fun_RLogEMean(Train,10);
        R_LogCC = fun_LogCC_new(Train,R_KA,10);
        %����ź�
        x0=alpha(m)*s+x0;%+pp;    %%%%%%%  ��Ҫ  %%%%%%%%%%%%%
        %�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ANMF_LogM
        T_LogM = fun_ANMF(R_LogM,x0,s);
        if T_LogM>1
            T_LogM = 1;
        end
        %%%% ANMF_LogCC
        T_LogCC = fun_ANMF(R_LogCC,x0,s);
        %%�ж�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if T_LogM>Th_LogM;          counter_LogM=counter_LogM+1;    end
        if T_LogCC>Th_LogCC;        counter_LogCC=counter_LogCC+1;    end
    end
    Pd_LogM_mc(m)=counter_LogM/MonteCarloPd;        counter_LogM=0;
    Pd_LogCC_mc(m)=counter_LogCC/MonteCarloPd;      counter_LogCC=0;
end
toc
close(h)
figure();
hold on

plot(SNRout,Pd_LogM_mc,'r.-','linewidth',1)
plot(SNRout,Pd_LogCC_mc,'r-*','linewidth',1)

h_leg = legend('ANMF with LogM','ANMF with LogCC');

xlabel('SNR/dB','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',20)
% set(h_leg,'Location','SouthEast')
grid on
% grid minor
str = ['PDLogM_',num2str(L),'Second','_s',num2str(sigma_t),'_',str_train,'.mat'];
save (str); 
