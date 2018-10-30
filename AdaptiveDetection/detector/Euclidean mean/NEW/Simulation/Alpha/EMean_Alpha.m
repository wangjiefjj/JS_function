clc
clear 
close all
% warning off
n = 2; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
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
    alpha_T = zeros(1,L_R);
    alpha_S = zeros(1,L_R);
    alpha_P = zeros(1,L_R);
    alpha_ML = zeros(1,L_R);
    parfor i =1:L_R
        R_KA = zeros(size(R));
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
        R_KA = R.*(t*t');
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
        R_KA2 = tau0^2*R.*(t*t');
        R_SCM = (fun_SCMN(Train));
        [R_CC,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
        [R_ECCT,alpha_T(i)] = fun_PowerCC(Train,R_KA,1,4);
        [R_ECCS,alpha_S(i)] = fun_PowerCC(Train,R_KA,1,8);
        [R_ECCP,alpha_P(i)] = fun_PowerCC(Train,R_KA,1,7);
        [R_ML,alpha_ML(i)]= fun_MLalpha(Train,R_SCM,R_KA,x0);
    end
    m_alpha(i_s) = mean(alpha);
    m_alpha_T(i_s) = mean(alpha_T);
    m_alpha_S(i_s) = mean(alpha_S);
    m_alpha_P(i_s) = mean(alpha_P);
    m_alpha_ML(i_s) = mean(alpha_ML);
end
close(h)
figure
hold on 
plot(sigma_t,m_alpha,'b','LineWidth',2)%k-.
plot(sigma_t,m_alpha_ML,'c','LineWidth',2)%k:
plot(sigma_t,m_alpha_T,'k','LineWidth',2)
plot(sigma_t,m_alpha_S,'k-*','LineWidth',2)
plot(sigma_t,m_alpha_P,'k-o','LineWidth',2)
h_leg = legend('CC','ML','ECCT','ECCS','ECCP');
grid on
box on
xlabel('\sigma^2','FontSize',10)
ylabel('\alpha','FontSize',10)
set(gca,'FontSize',10)
set(h_leg,'Location','SouthEast')
% str = [str_train,'_','Alpha','_',num2str(n),'N','.mat'];
% save (str); 


