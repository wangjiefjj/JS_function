clc
clear 
close all
% warning off
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�������ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG��������ͬ
sigma_t = 0.1;
rou = 0.95;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
R = fun_rho(rou,N,2);
SNRout=-5:1:20; % ���SNR
PFA=1e-3;% PFA=1e-4;
SNRnum=10.^(SNRout/10);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.1;
nn = 0:N-1;
s = exp(-1i*2*pi*nn*theta_sig)'; %%%%%% ϵͳ����ʸ��
R_KA = zeros(size(R));
for i = 1:1000
    t = normrnd(1,sigma_t,N,1);%%0~0.5%%ʧ������
    R_KA = R_KA + R.*(t*t')/1000;
end
iter = 10;
for i =1:1000
    warning off
    i
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_x0 = (fun_SCMN(x0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [R_result,a(i)] = fun_LogCC(Train,R_SCM,R_KA);
     R_LogMean = fun_RLogEMean(Train,2);
    [R_LCC,alcc(i)] = fun_test_LogEKA(Train,R_SCM,R_KA);
    [R_result,a(i)] = fun_CCIter(Train,R_LogMean,R_KA);
    [R_CC,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
%     R_2 = 0.5 * R_KA + 0.5 * R_SCM;
    error_RLCC(i) = norm(R_LCC-R,'fro')/norm(R,'fro');
    error_RL(i) = norm(R_LogMean-R,'fro')/norm(R,'fro');
    error_R(i) = norm(R_result-R,'fro')/norm(R,'fro');
    error_RCC(i) = norm(R_CC-R,'fro')/norm(R,'fro');
%     error_R_2(i) = norm(R_2-R,'fro')/norm(R,'fro');
    error_RSCM(i) = norm(R_SCM-R,'fro')/norm(R,'fro');
end
m_errorRLCC = mean(error_RLCC)/norm(R,'fro');
m_errorRL = mean(error_RL)/norm(R,'fro');
m_errorR = mean(error_R)/norm(R,'fro');
m_errorRCC = mean(error_RCC)/norm(R,'fro');
% m_errorR_2 = mean(error_R_2)/norm(R,'fro');
m_errorRSCM = mean(error_RSCM)/norm(R,'fro');
m_errorRKA = norm(R_KA-R,'fro')/norm(R,'fro');
mean_alcc = mean(alcc);
mean_a = mean(a);
mean_alpha = mean(alpha);
% plot(a,'r')
% hold on
% plot(alpha,'k')
% axis([1,1000,0,1]);
