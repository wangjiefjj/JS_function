clc
clear 
close all
% warning off
n = 1; %几倍的样本
str_train = 'g';%%训练数据分布，p:IG纹理复合高斯，k：k分布，g：gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG的选项，1为每个距离单元IG纹理都不同
sigma_t = [0.01,0.1:0.1:10];
%  sigma_t = [0.01:0.01:1];
% sigma_t = [11:101];
L_s = length(sigma_t);
L_R = 1000;
opt = 'k';
rou = 0.95;  %%协方差矩阵生成的迟滞因子
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
R = fun_rho(rou,N,2);
L=round(n*N); 
h = waitbar(0,'Please wait...');
for i_s = 1:L_s
    waitbar(i_s/L_s,h,sprintf([num2str(i_s/L_s*100),'%%']));
    R_KA = zeros(size(R));
    for i = 1:1000
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%失配向量
        R_KA = R_KA + R.*(t*t')/1000;
    end
    parfor i =1:L_R
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
        R_SCM = (fun_SCMN(Train));
        R_NSCM = (fun_NSCMN(Train));
        R_x0 = (fun_SCMN(x0));
        R_AML = fun_AML(Train);
        [R_CC,alpha(i)]=fun_CC(Train,R_SCM,R_KA);
        [R_CCIter,alpha_iter(i)]=fun_CCIter(Train,R_SCM,R_KA);
%         if sigma_t(i_s) < 0.2
%             [R_CCIter,alpha_iter(i)]=fun_CCIter2(Train,R_SCM,R_KA);
%         else
%             [R_CCIter,alpha_iter(i)]=fun_CCIter(Train,R_SCM,R_KA);
%         end
%         if sigma_t(i_s)<0.9
%            [R_AMLCC,alpha_aml(i)]=fun_AMLCC5(Train,R_KA);
%         else
%            [R_AMLCC,alpha_aml(i)]=fun_CCIter(Train,R_AML,R_KA);
%         end
        [R_CCML,alpha_ML(i)]=fun_MLalpha(Train,R_SCM,R_KA,x0);
    end
    m_alpha(i_s) = mean(alpha);
    m_alpha_iter(i_s) = mean(alpha_iter);
    m_alpha_ML(i_s) = mean(alpha_ML);
end
close(h)
str = [str_train,'_','CCIter_Alpha','_',num2str(n),'N','.mat'];
save (str); 


