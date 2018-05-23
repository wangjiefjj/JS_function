clc
clear 
close all
warning off
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = 0.1;
opt = 'k';
rou = 0.95;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
Sigma = fun_rho(rou,N,2);
SNRout=-5:1:20; % ���SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
R_KA = zeros(size(Sigma));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
    R_KA = R_KA + Sigma.*(t*t')/1000;
end
tic
for i =1:1000
    i
    Train = fun_TrainData(str_train,N,L,Sigma,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    [x0,tau] = fun_TrainData(str_train,N,1,Sigma,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [R_CC,alpha_cc(i)]=fun_CC(Train,R_SCM,R_KA);
    [R_ML,alpha_ML(i)]=fun_MLalpha(Train,R_SCM,R_KA,x0);
    [R_ECC,alpha_ecc(i)]=fun_ECC(Train,R_SCM,R_KA);
    error_RCC(i) = norm(R_CC-Sigma,'fro');
    error_RECC(i) = norm(R_ECC-Sigma,'fro');
    error_RML(i) = norm(R_ML-Sigma,'fro');
    error_RSCM(i) = norm(R_SCM-Sigma,'fro');
    error_RNSCM(i) = norm(R_NSCM-Sigma,'fro');
end
toc
m_errorRCC = mean(error_RCC)/norm(Sigma,'fro');
m_errorRECC = mean(error_RECC)/norm(Sigma,'fro');
m_errorRML = mean(error_RML)/norm(Sigma,'fro');
m_errorRSCM = mean(error_RSCM)/norm(Sigma,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(Sigma,'fro');
m_errorRKA = norm(R_KA-Sigma,'fro')/norm(Sigma,'fro');

mean_alpha_ecc = mean(alpha_ecc);
mean_alpha_cc = mean(alpha_cc);
mean_alpha_ML = mean(alpha_ML);

num = 1:1000;
plot(num,alpha_cc)
hold on
plot(num,alpha_ecc)
plot(num,alpha_ML)
axis([1,1000,0,2])




