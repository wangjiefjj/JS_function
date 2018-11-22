%%选择数据的最佳模型
%%SCNR,K,固定，SIRP环境，精度随参数lambda的改变
clc
clear
close all
Class=3; %%
rho=4;  %%GIC的参数 
MC = 1000;
rou = 0.90;  %%协方差矩阵生成的迟滞因子
fc = 0;
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
lambda = 1;%%%越小非高斯越严重
mu = 1:0.1:3;
% n =1.1:0.3:8; %几倍的样本
% L=round(n*N);SNRout=10;
SCNRout=0;
SCNRnum=10.^(SCNRout/10);
CNRout=30;
CNRnum=10.^(CNRout/10); 
L = 20;
theta_sig = 0.2;
nn = 0:N-1;
p = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
str_train = 'p';
%%%IG的选项，1为每个距离单元IG纹理都不同
%2每个单元纹理值一样，为部分均匀环境
opt_train = 1;     
%%杂波协方差
Rc1 = fun_rho(rou,N,1,fc);
Rc1 = CNRnum * Rc1;
Rc1 = Rc1+ eye(N) ;%+ eye(N)
Rc2 = Rc1;
str=['SIRP_mu_L',num2str(L),'_rho',num2str(rho),'.mat'];
 a = sqrt(SCNRnum/abs(p'/Rc2*p));
tic
h = waitbar(1,'Please wait...');
for i_mu = 1:length(mu)   
    waitbar(i_mu/length(mu),h,sprintf([num2str(i_mu/length(mu)*100),'%%']));
    count_AIC1 = 0;
    count_ABIC1 = 0;
    count_GIC1 = 0;
    count_AICc1 = 0;


    count_AIC2 = 0;
    count_ABIC2 = 0;
    count_GIC2 = 0;
    count_AICc2 = 0;
    
    count_AIC3 = 0;
    count_ABIC3 = 0;
    count_GIC3 = 0;
    count_AICc3 = 0;
    %%参数个数
    %只用辅助数据时
    H1_num1 = N^2;
    H2_num1 = N^2+1;
    H3_num1 = N^2+L;
    %主辅数据时
    H1_num2 = N^2+1;
    H2_num2 = N^2+3;
    H3_num2 = N^2+L+2;
    %主数据时
    H1_num3 = N^2+1;
    H2_num3 = N^2+2;
    H3_num3 = N^2+2;
    parfor i = 1:MC
        warning off
        %%产生数据
        [Train,tauk] = fun_TrainData(str_train,N,L,Rc1,lambda,mu(i_mu),opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData('g',N,1,Rc2,lambda,mu(i_mu),opt_train); % 接收信号仅包括杂波和噪声
        x0 =x0+a*p;
        %%s函数计算 
        %只用辅助数据时
        s_H1_1 = -abs(fun_s_H1(Train,p,1));
        s_H2_1 = -abs(fun_s_H2(Train,p,1));
        [s_H3_1] = -abs(fun_s_H3(Train,p,1));
        %主辅数据时
        s_H1_2 = -abs(fun_s_H1([Train,x0],p,2));
        s_H2_2 = -abs(fun_s_H2([Train,x0],p,2));
        s_H3_2 = -abs(fun_s_H3([Train,x0],p,2));
        %只用主数据时
        s_H1_3 = -abs(fun_s_H1([Train,x0],p,3));
        s_H2_3 = -abs(fun_s_H2([Train,x0],p,3));
        s_H3_3 = -abs(fun_s_H3([Train,x0],p,3));
        
        %%模型选择准则%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %只用辅助数据时%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%AIC%%%%%%%%%
        H1_AIC_1 = fun_AIC(s_H1_1,H1_num1);
        H2_AIC_1 = fun_AIC(s_H2_1,H2_num1);
        H3_AIC_1 = fun_AIC(s_H3_1,H3_num1);
        Class_AIC1 =[H1_AIC_1,H2_AIC_1,H3_AIC_1];
        [~,Class_AIC_num1] = min(Class_AIC1);
        if Class_AIC_num1 == Class  
            count_AIC1 = count_AIC1+1;  
        end
        %%GIC%%%%%%
        H1_GIC_1 = fun_GIC(s_H1_1,H1_num1,rho);
        H2_GIC_1 = fun_GIC(s_H2_1,H2_num1,rho);
        H3_GIC_1 = fun_GIC(s_H3_1,H3_num1,rho);
        Class_GIC1 =[H1_GIC_1,H2_GIC_1,H3_GIC_1];
        [~,Class_GIC_num1] = min(Class_GIC1);
        if Class_GIC_num1 == Class  
            count_GIC1 = count_GIC1+1;  
        end
        %%AICc%%%%%%%
        H1_AICc_1 = fun_AICc(s_H1_1,H1_num1,N,L);
        H2_AICc_1 = fun_AICc(s_H2_1,H2_num1,N,L);
        H3_AICc_1 = fun_AICc(s_H3_1,H3_num1,N,L);
        Class_AICc1 =[H1_AICc_1,H2_AICc_1,H3_AICc_1];
        [~,Class_AICc_num1] = min(Class_AICc1);
        if Class_AICc_num1 == Class  
            count_AICc1 = count_AICc1+1;  
        end
        %%ABIC%%%%%%%
        H1_ABIC_1 = fun_ABIC(s_H1_1,H1_num1,L);
        H2_ABIC_1 = fun_ABIC(s_H2_1,H2_num1,L);
        H3_ABIC_1 = fun_ABIC(s_H3_1,H3_num1,L);
        Class_ABIC1 =[H1_ABIC_1,H2_ABIC_1,H3_ABIC_1];
        [~,Class_ABIC_num1] = min(Class_ABIC1);
        if Class_ABIC_num1 == Class  
            count_ABIC1 = count_ABIC1+1;  
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %主辅数据时%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%AIC%%%%%%%%%%%%%%%%
        AIC_H1_2 = fun_AIC(s_H1_2,H1_num2);
        AIC_H2_2 = fun_AIC(s_H2_2,H2_num2);
        AIC_H3_2 = fun_AIC(s_H3_2,H3_num2);
        Class_AIC2 =[AIC_H1_2,AIC_H2_2,AIC_H3_2];
        [~,Class_AIC_num2] = min(Class_AIC2);
        if Class_AIC_num2 == Class  
            count_AIC2 = count_AIC2+1;  
        end
        %%GIC%%%%%%%%%%%%%%
        H1_GIC_2 = fun_GIC(s_H1_2,H1_num2,rho);
        H2_GIC_2 = fun_GIC(s_H2_2,H2_num2,rho);
        H3_GIC_2 = fun_GIC(s_H3_2,H3_num2,rho);
        Class_GIC2 =[H1_GIC_2,H2_GIC_2,H3_GIC_2];
        [~,Class_GIC_num2] = min(Class_GIC2);
        if Class_GIC_num2 == Class  
            count_GIC2 = count_GIC2+1;  
        end
        %%AICc%%%%%%%
        H1_AICc_2 = fun_AICc(s_H1_2,H1_num2,N,L);
        H2_AICc_2 = fun_AICc(s_H2_2,H2_num2,N,L);
        H3_AICc_2 = fun_AICc(s_H3_2,H3_num2,N,L);
        Class_AICc2 =[H1_AICc_2,H2_AICc_2,H3_AICc_2];
        [~,Class_AICc_num2] = min(Class_AICc2);
        if Class_AICc_num2 == Class  
            count_AICc2 = count_AICc2+1;  
        end
        %%ABIC%%%%%%%
        H1_ABIC_2 = fun_ABIC(s_H1_2,H1_num2,L);
        H2_ABIC_2 = fun_ABIC(s_H2_2,H2_num2,L);
        H3_ABIC_2 = fun_ABIC(s_H3_2,H3_num2,L);
        Class_ABIC2 =[H1_ABIC_2,H2_ABIC_2,H3_ABIC_2];
        [~,Class_ABIC_num2] = min(Class_ABIC2);
        if Class_ABIC_num2 == Class  
            count_ABIC2 = count_ABIC2+1;  
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %主数据时%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%AIC%%%%%%%%%%%%%%%%
        AIC_H1_3 = fun_AIC(s_H1_3,H1_num3);
        AIC_H2_3 = fun_AIC(s_H2_3,H2_num3);
        AIC_H3_3 = fun_AIC(s_H3_3,H3_num3);
        Class_AIC3 =[AIC_H1_3,AIC_H2_3,AIC_H3_3];
        [~,Class_AIC_num3] = min(Class_AIC3);
        if Class_AIC_num3 == Class  
            count_AIC3 = count_AIC3+1;  
        end
        %%GIC%%%%%%%%%%%%%%
        H1_GIC_3 = fun_GIC(s_H1_3,H1_num3,rho);
        H2_GIC_3 = fun_GIC(s_H2_3,H2_num3,rho);
        H3_GIC_3 = fun_GIC(s_H3_3,H3_num3,rho);
        Class_GIC3 =[H1_GIC_3,H2_GIC_3,H3_GIC_3];
        [~,Class_GIC_num3] = min(Class_GIC3);
        if Class_GIC_num3 == Class  
            count_GIC3 = count_GIC3+1;  
        end
        %%AICc%%%%%%%
        H1_AICc_3 = fun_AICc(s_H1_3,H1_num3,N,L);
        H2_AICc_3 = fun_AICc(s_H2_3,H2_num3,N,L);
        H3_AICc_3 = fun_AICc(s_H3_3,H3_num3,N,L);
        Class_AICc3 =[H1_AICc_3,H2_AICc_3,H3_AICc_3];
        [~,Class_AICc_num3] = min(Class_AICc3);
        if Class_AICc_num3 == Class  
            count_AICc3 = count_AICc3+1;  
        end
        %%ABIC%%%%%%%
        H1_ABIC_3 = fun_ABIC(s_H1_3,H1_num3,L);
        H2_ABIC_3 = fun_ABIC(s_H2_3,H2_num3,L);
        H3_ABIC_3 = fun_ABIC(s_H3_3,H3_num3,L);
        Class_ABIC3 =[H1_ABIC_3,H2_ABIC_3,H3_ABIC_3];
        [~,Class_ABIC_num3] = min(Class_ABIC3);
        if Class_ABIC_num3 == Class  
            count_ABIC3 = count_ABIC3+1;  
        end
    end
    Accuracy_AIC1(i_mu) = count_AIC1/MC;
    Accuracy_GIC1(i_mu) = count_GIC1/MC;
    Accuracy_AICc1(i_mu) = count_AICc1/MC;
    Accuracy_ABIC1(i_mu) = count_ABIC1/MC;

    Accuracy_AIC2(i_mu) = count_AIC2/MC;
    Accuracy_GIC2(i_mu) = count_GIC2/MC;
    Accuracy_AICc2(i_mu) = count_AICc2/MC;
    Accuracy_ABIC2(i_mu) = count_ABIC2/MC;
    
    Accuracy_AIC3(i_mu) = count_AIC3/MC;
    Accuracy_GIC3(i_mu) = count_GIC3/MC;
    Accuracy_AICc3(i_mu) = count_AICc3/MC;
    Accuracy_ABIC3(i_mu) = count_ABIC3/MC;
end
toc
close(h)
% figure(1)
% hold on
% plot(SNRout,Accuracy_AIC1,'r-s','LineWidth',2)
% plot(SNRout,Accuracy_GIC1,'g-o','LineWidth',2)
% plot(SNRout,Accuracy_AICc1,'b-*','LineWidth',2)
% plot(SNRout,Accuracy_ABIC1,'k-.','LineWidth',2)
% title('只用辅助数据')
% h_leg1 = legend('AIC','GIC','AICc','ABIC');
% xlabel('\mu')
% ylabel('Accuracy')
% set(gca,'FontSize',10)
% set(h_leg1,'Location','SouthEast')
% grid on
% box on
figure(2)
hold on
plot(mu,Accuracy_AIC2,'r-s','LineWidth',2)
plot(mu,Accuracy_GIC2,'g-o','LineWidth',2)
plot(mu,Accuracy_AICc2,'b-*','LineWidth',2)
plot(mu,Accuracy_ABIC2,'k-.','LineWidth',2)
title('主辅数据')
h_leg2 = legend('AIC','GIC','AICc','ABIC');
xlabel('\mu')
ylabel('Accuracy')
set(gca,'FontSize',10)
set(h_leg2,'Location','SouthEast')
grid on
box on
figure(3)
hold on
plot(mu,Accuracy_AIC3,'r-s','LineWidth',2)
plot(mu,Accuracy_GIC3,'g-o','LineWidth',2)
plot(mu,Accuracy_AICc3,'b-*','LineWidth',2)
plot(mu,Accuracy_ABIC3,'k-.','LineWidth',2)
title('主数据')
h_leg2 = legend('AIC','GIC','AICc','ABIC');
xlabel('\mu')
ylabel('Accuracy')
set(gca,'FontSize',10)
set(h_leg2,'Location','SouthEast')
grid on
box on
save(str)
