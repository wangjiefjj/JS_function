clc
clear 
close all
% warning off
n = 0.5; %����������
str_train = 'g';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 2;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = 0.1;
rou = 0.95;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
CNRout=10; % ���CNR
PFA=1e-3;% PFA=1e-4;
CNRnum=10.^(CNRout/10);
R = fun_rho(rou,N,1,0.1);
R = CNRnum*R + eye(N);
MonteCarloPfa=1/PFA*100;
MonteCarloPd=1e4;
L=round(n*N); 
theta_sig = 0.5;
nn = 0:N-1;

iter = 10;
for i =1:1000
    warning off
    i;
    Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
    x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
    R_SCM = (fun_SCMN(Train));
    R_SCM2 = R_SCM + eye(N);
    R_NSCM = (fun_NSCMN(Train));
    R_x0 = (fun_SCMN(x0));
    R_MESS = fun_MESS(R_SCM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_RSCM(i) = norm(R_SCM-R,'fro');
    error_RSCM2(i) = norm(R_SCM2-R,'fro');
    error_RNSCM(i) = norm(R_NSCM-R,'fro');
    error_RMESS(i) = norm(R_MESS-R,'fro');   
end
m_errorRSCM = mean(error_RSCM)/norm(R,'fro');
m_errorRSCM2 = mean(error_RSCM2)/norm(R,'fro');
m_errorRNSCM = mean(error_RNSCM)/norm(R,'fro');
m_errorRMESS = mean(error_RMESS)/norm(R,'fro');
% rank(R_SCM)
% rank(R_MESS)