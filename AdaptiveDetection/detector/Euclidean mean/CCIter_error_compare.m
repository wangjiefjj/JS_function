clc
clear 
close all
% warning off
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = sqrt([0.01,0.1:0.1:0.9]);
L_s = length(sigma_t);
L_R = 100;
opt = 'k';
rou = 0.95;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
R = fun_rho(rou,N,2);
L=round(n*N); 
m_errorRCC = zeros(1,L_s);
m_errorRCCIter = zeros(1,L_s);
m_errorRCCML= zeros(1,L_s);
m_errorR_2 = zeros(1,L_s);
m_errorRSCM = zeros(1,L_s);
m_errorRKA = zeros(1,L_s);
for i_s = 1:L_s
    i_s
    error_RCC = zeros(1,L_R);
    error_RECCT = zeros(1,L_R);
    error_RML = zeros(1,L_R);
    error_R_2 = zeros(1,L_R);
    error_RSCM = zeros(1,L_R);
    
    parfor i =1:L_R
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        [x0,tau0] = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
%         R_real = (tau0^2)*R;
        %%����Э����
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%ʧ������
        R_KA =  (tau0^2*R).*(t*t');
        %%%
        %%%%Э�������%%%%%%%%%%%%%%%%%%%%%%    
        R_SCM = (fun_SCMN(Train));
        R_CC = fun_CC(Train,R_SCM,R_KA);
        R_ECCT = fun_PowerCC(Train,R_KA,1,4);
        R_ECCS = fun_PowerCC(Train,R_KA,1,8);
        R_ECCP = fun_PowerCC(Train,R_KA,1,7);
        R_ML = fun_MLalpha(Train,R_SCM,R_KA,x0);
        %%NFN
        error_RCC(i) = norm(R_CC-R_real,'fro')/norm(R,'fro');
        error_RECCT(i) = norm(R_ECCT-R_real,'fro')/norm(R,'fro');
        error_RECCS(i) = norm(R_ECCS-R_real,'fro')/norm(R,'fro');
        error_RECCP(i) = norm(R_ECCP-R_real,'fro')/norm(R,'fro');
        error_RML(i) = norm(R_ML-R_real,'fro')/norm(R_real,'fro');    
        error_RSCM(i) = norm(R_SCM-R_real,'fro')/norm(R,'fro');
    end

    m_errorRCC(i_s) = mean(error_RCC);
    m_errorRCCIter(i_s) = mean(error_RECCT);
%     m_errorRLECC(i_s) = mean(error_RLECC);
    m_errorRCCML(i_s) = mean(error_RML);
    m_errorRSCM(i_s) = mean(error_RSCM);
    

    m_alpha(i_s) = mean(alpha);
    m_alpha_iter(i_s) = mean(alpha_iter);
%     m_alpha_lecc(i_s) = mean(alpha_lecc);
    m_alpha_ML(i_s) = mean(alpha_ML);
end

figure
hold on 
plot(sigma_t,m_errorRSCM,'r','LineWidth',2)%k--
plot(sigma_t,m_errorRCC,'b','LineWidth',2)%k-.
plot(sigma_t,m_errorRCCIter, 'k','LineWidth',2)
plot(sigma_t,m_errorRCCML,'g','LineWidth',2)%k:
% plot(sigma_t,m_errorRLECC,'y','LineWidth',2)%k:
grid on
h_leg = legend('SCM','CC','KA-CE','ML');
xlabel('\sigma^2','FontSize',10)
ylabel('Error','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')
% axis([0,10,0,100])
