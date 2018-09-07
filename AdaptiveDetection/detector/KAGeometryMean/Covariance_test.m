clc
clear 
close all
warning off
n = 0.5; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = 0.9;
rou = 0.90;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
Sigma = fun_rho(rou,N,1);
SNRout=-5:1:20; % ���SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'/sqrt(N); %%%%%% ϵͳ����ʸ��
% R_KA = zeros(size(Sigma));
% for i = 1:1000
%     t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
%     R_KA = R_KA + Sigma.*(t*t')/1000;
% end
tic
RR = zeros(N,N);
parfor i =1:1e3
    Train = fun_TrainData(str_train,N,L,Sigma,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    [x0,tau] = fun_TrainData(str_train,N,1,Sigma,lambda,mu,opt_train); 
    R_KA = zeros(size(Sigma));
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
    R_KA = Sigma.*(t*t');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_LogM = fun_RLogEMean(Train,4);
%     [R_CC,alpha_cc(i)]=fun_CC(Train,R_NSCM,R_KA);
%     [R_ML,alpha_ML(i)]=fun_MLalpha(Train,R_NSCM,R_KA,x0);
%     [R_ECC,alpha_ecc(i)]=fun_PowerCC(Train,R_KA,1,4);
    [R_LogCC,alpha_lecc(i)]=fun_LogCC_new(Train,R_KA,4);
    [R_ECC,alpha_ecc(i)]=fun_PowerCC(Train,R_KA,1,4);
    [R_PCC,alpha_pcc(i)]=fun_PowerCC(Train,R_KA,-1,4);
    R_SFP = fun_SFP(Train,1);
    error_PCC(i) = norm(R_PCC-Sigma,'fro');
    error_LogCC(i) = norm(R_LogCC-Sigma,'fro');
    error_LogM(i) = norm((R_LogM)-(Sigma),'fro');
    error_RSFP(i) = norm(R_SFP-Sigma,'fro');
    ANMF_LogM(i) = (fun_ANMF(R_LogM,x0,s));
%     ANMF_LogM(i)
%     if ANMF_LogM(i)>1
%         ANMF_LogM(i)
%         RR = R_LogM;
%         break
%     end
end
toc
m_errorLogCC = mean(error_LogCC)/norm(Sigma,'fro');
m_errorPCC = mean(error_PCC)/norm(Sigma,'fro');
m_errorLogM = mean(error_LogM)/norm(Sigma,'fro');
m_errorRSFP = mean(error_RSFP)/norm(Sigma,'fro');
% 
% 
% mean_alpha_ecc = mean(alpha_ecc);
mean_alpha_lecc = mean(alpha_lecc);
mean_alpha_ecc = mean(alpha_ecc);
mean_alpha_pcc = mean(alpha_pcc);
% mean_alpha_cc = mean(alpha_cc);
% mean_alpha_ML = mean(alpha_ML);
TANMF_LogM=sort(ANMF_LogM,'descend');

% num = 1:1000;
% plot(num,alpha_cc,'b')
% hold on
% plot(num,alpha_ML,'g')
% plot(num,alpha_ecc,'r')
% plot(num,alpha_lecc,'k')
% legend('CC','ML','E','LE')




