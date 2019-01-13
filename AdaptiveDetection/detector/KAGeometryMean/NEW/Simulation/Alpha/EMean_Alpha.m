clc
clear 
close all
% warning off
n = 0.5; %几倍的样本
str_train = 'p';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = sqrt([0.01,0.1:0.1:1]);
L_s = length(sigma_t);
L_R = 1000;
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
R = fun_rho(rou,N,1,0.05);
L=round(n*N); 
h = waitbar(0,'Please wait...');
for i_s = 1:L_s
    waitbar(i_s/L_s,h,sprintf([num2str(i_s/L_s*100),'%%']));
    alpha = zeros(1,L_R);
    alpha_ECC = zeros(1,L_R);
    alpha_PCC = zeros(1,L_R);
    alpha_LogCC = zeros(1,L_R);
    alpha_SFP = zeros(1,L_R);
    parfor i =1:L_R
        R_KA = zeros(size(R));
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
        R_KA = R.*(t*t');
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
        R_KA2 = tau0^2*R.*(t*t');
        R_SCM = (fun_SCMN(Train));
        [~,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
        [~,alpha_ECC(i)] = fun_PowerCC(Train,R_KA,1,10);
        [~,alpha_PCC(i)]=fun_PowerCC(Train,R_KA,-1,10);
        [~,alpha_LogCC(i)]=fun_LogCC_new(Train,R_KA,10);
        [~,alpha_SFP(i)] = fun_SFP(Train,1);
        
    end
    m_alpha_ECC(i_s) = mean(alpha_ECC);
    m_alpha_PCC(i_s) = mean(alpha_PCC);
    m_alpha_LogCC(i_s) = mean(alpha_LogCC);
    m_alpha_CC(i_s) = mean(alpha);
    m_alpha_SFP(i_s) = mean(alpha_SFP);
end
close(h)
figure
hold on 
plot(sigma_t.^2,m_alpha_ECC,'b','LineWidth',2)%k-.
plot(sigma_t.^2,m_alpha_PCC,'c','LineWidth',2)%k:
plot(sigma_t.^2,m_alpha_LogCC,'k','LineWidth',2)
plot(sigma_t.^2,m_alpha_CC,'y','LineWidth',2)
plot(sigma_t.^2,m_alpha_SFP,'r','LineWidth',2)
h_leg = legend('KA-E','KA-PE','KA-LogE','CC','SFP');
grid on
box on
xlabel('\sigma^2','FontSize',10)
ylabel('\alpha','FontSize',10)
axis([0.01,1,0,1])
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')
% str = [str_train,'_','Alpha','_',num2str(n),'N','.mat'];
% save (str); 


