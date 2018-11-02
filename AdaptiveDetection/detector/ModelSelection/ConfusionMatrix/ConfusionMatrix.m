%%计算混淆矩阵
clc
clear
close all
rho=2;  %%GIC的参数 
MC = 5000;
rou = 0.90;  %%协方差矩阵生成的迟滞因子
fc = 0;
%%%%假设参数设置
Na = 2;     % 阵元数
Np = 4;     % 脉冲数
N = Na*Np;
lambda = 1;%%%越小非高斯越严重
mu = 1;
SNRout=10;
SNRnum=10.^(SNRout/10);
CNRout=30;
CNRnum=10.^(CNRout/10);
L = round(1.1*N);
theta_sig = 0.2;
nn = 0:N-1;
p = exp(-1i*2*pi*nn*theta_sig).'/sqrt(N); %%%%%% 系统导向矢量
a = sqrt(SNRnum);

for Class = 1:3
    if Class==1%%均匀
        str_train = 'g';
        %%杂波协方差
        sigmaf = 0.03; %%杂波谱展宽
        rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
        Rc1 = CNRnum * toeplitz(rc);
        Rc1 = Rc1+ eye(N) ;%
        Rc2 = Rc1;
        opt_train = 0;
    elseif Class == 2%%部分均匀
        str_train = 'p';
        opt_train = 2;
        %%杂波协方差
        sigmaf = 0.03; %%杂波谱展宽
        rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
        Rc1 = CNRnum * toeplitz(rc);
        Rc1 = Rc1+ eye(N) ;%+ eye(N)
        Rc2 = Rc1;
    elseif Class == 3%%SIRP
        str_train = 'p';
        %%%IG的选项，1为每个距离单元IG纹理都不同
        %2每个单元纹理值一样，为部分均匀环境
        opt_train = 1;     
        %%杂波协方差
        sigmaf = 0.03; %%杂波谱展宽
        rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
        Rc1 = CNRnum * toeplitz(rc);
        Rc1 = Rc1+ eye(N) ;%+ eye(N)
        Rc2 = Rc1;
    end
    %%
    %%三种环境计数
    count_AIC_H1 = 0;
    count_AIC_H2 = 0;
    count_AIC_H3 = 0;
    
    count_ABIC_H1 = 0;
    count_ABIC_H2 = 0;
    count_ABIC_H3 = 0;
    
    count_GIC1_H1 = 0;
    count_GIC1_H2 = 0;
    count_GIC1_H3 = 0;
    
    count_GIC2_H1 = 0;
    count_GIC2_H2 = 0;
    count_GIC2_H3 = 0;
    
    count_AICc_H1 = 0;
    count_AICc_H2 = 0;
    count_AICc_H3 = 0;
    %%
    %MC
    parfor i = 1:MC
        warning off
        %%产生数据
        [Train,tauk] = fun_TrainData(str_train,N,L,Rc1,lambda,mu,opt_train);%%产生的训练数据,协方差矩阵为rouR的高斯杂波
        [x0,tau0] = fun_TrainData(str_train,N,1,Rc2,lambda,mu,opt_train); % 接收信号仅包括杂波和噪声
        x0 =x0+a*p;
        %主辅数据时
        s_H1 = -abs(fun_s_H1([Train,x0],p,2));
        s_H2 = -abs(fun_s_H2([Train,x0],p,2));
        s_H3 = -abs(fun_s_H3([Train,x0],p,2));
        
        %%参数个数
        %主辅数据时
        H1_num = N^2+1;
        H2_num = N^2+3;
        H3_num = N^2+L+2;
        %%
        %%模型选择准则%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %主数据时%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        %%AIC%%%%%%%%%%%%%%%%
        AIC_H1 = fun_AIC(s_H1,H1_num);
        AIC_H2 = fun_AIC(s_H2,H2_num);
        AIC_H3 = fun_AIC(s_H3,H3_num);
        Class_AIC =[AIC_H1,AIC_H2,AIC_H3];
        [~,Class_AIC_num] = min(Class_AIC);
        if Class_AIC_num == 1  
            count_AIC_H1 = count_AIC_H1+1;  
        elseif Class_AIC_num == 2  
            count_AIC_H2 = count_AIC_H2+1;
        elseif Class_AIC_num == 3  
            count_AIC_H3 = count_AIC_H3+1;    
        end
       %%
        %%GIC rho=2%%%%%%%%%%%%%%
        H1_GIC1 = fun_GIC(s_H1,H1_num,2);
        H2_GIC1 = fun_GIC(s_H2,H2_num,2);
        H3_GIC1 = fun_GIC(s_H3,H3_num,2);
        Class_GIC1 =[H1_GIC1,H2_GIC1,H3_GIC1];
        [~,Class_GIC1_num1] = min(Class_GIC1);
        if Class_GIC1_num1 == 1  
            count_GIC1_H1 = count_GIC1_H1+1;  
        elseif Class_GIC1_num1 == 2  
            count_GIC1_H2 = count_GIC1_H2+1; 
        elseif Class_GIC1_num1 == 3  
            count_GIC1_H3 = count_GIC1_H3+1;
        end
        %%
        %%GIC rho=4%%%%%%%%%%%%%%
        H1_GIC2 = fun_GIC(s_H1,H1_num,4);
        H2_GIC2 = fun_GIC(s_H2,H2_num,4);
        H3_GIC2 = fun_GIC(s_H3,H3_num,4);
        Class_GIC2 =[H1_GIC2,H2_GIC2,H3_GIC2];
        [~,Class_GIC2_num1] = min(Class_GIC2);
        if Class_GIC2_num1 == 1  
            count_GIC2_H1 = count_GIC2_H1+1;  
        elseif Class_GIC2_num1 == 2  
            count_GIC2_H2 = count_GIC2_H2+1; 
        elseif Class_GIC2_num1 == 3  
            count_GIC2_H3 = count_GIC2_H3+1;
        end
        %%
        %%AICc%%%%%%%
        H1_AICc = fun_AICc(s_H1,H1_num,N,L);
        H2_AICc = fun_AICc(s_H2,H2_num,N,L);
        H3_AICc = fun_AICc(s_H3,H3_num,N,L);
        Class_AICc =[H1_AICc,H2_AICc,H3_AICc];
        [~,Class_AICc_num] = min(Class_AICc);
        if Class_AICc_num == 1
            count_AICc_H1 = count_AICc_H1+1; 
        elseif Class_AICc_num == 2
            count_AICc_H2 = count_AICc_H2+1;
        elseif Class_AICc_num == 3    
            count_AICc_H3 = count_AICc_H3+1;  
        end
        %%
        %%ABIC%%%%%%%
        H1_ABIC = fun_ABIC(s_H1,H1_num,L);
        H2_ABIC = fun_ABIC(s_H2,H2_num,L);
        H3_ABIC = fun_ABIC(s_H3,H3_num,L);
        Class_ABIC =[H1_ABIC,H2_ABIC,H3_ABIC];
        [~,Class_ABIC_num] = min(Class_ABIC);
        if Class_ABIC_num == 1
            count_ABIC_H1 = count_ABIC_H1+1;
        elseif Class_ABIC_num == 2
            count_ABIC_H2 = count_ABIC_H2+1;
        elseif Class_ABIC_num == 3    
            count_ABIC_H3 = count_ABIC_H3+1;  
        end
    end
    %%
    %%Confusion Matrix compute
    Confusion_AIC(Class,1) = count_AIC_H1/MC;
    Confusion_AIC(Class,2) = count_AIC_H2/MC;
    Confusion_AIC(Class,3) = count_AIC_H3/MC;
    
    Confusion_GIC1(Class,1) = count_GIC1_H1/MC;
    Confusion_GIC1(Class,2) = count_GIC1_H2/MC;
    Confusion_GIC1(Class,3) = count_GIC1_H3/MC;
    
    Confusion_GIC2(Class,1) = count_GIC2_H1/MC;
    Confusion_GIC2(Class,2) = count_GIC2_H2/MC;
    Confusion_GIC2(Class,3) = count_GIC2_H3/MC;
    
    Confusion_AICc(Class,1) = count_AICc_H1/MC;
    Confusion_AICc(Class,2) = count_AICc_H2/MC;
    Confusion_AICc(Class,3) = count_AICc_H3/MC;
    
    Confusion_ABIC(Class,1) = count_ABIC_H1/MC;
    Confusion_ABIC(Class,2) = count_ABIC_H2/MC;
    Confusion_ABIC(Class,3) = count_ABIC_H3/MC;
end
   
str=['Confusion','_L',num2str(L),'.mat'];
save(str)
