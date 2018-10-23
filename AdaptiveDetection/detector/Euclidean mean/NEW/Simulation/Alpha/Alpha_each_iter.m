clc
clear 
close all
% warning off
n = 1; %����������
str_train = 'p';%%ѵ�����ݷֲ���p:IG�����ϸ�˹��k��k�ֲ���g��gauss
lambda = 3;
mu = 1;
opt_train = 1; %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
sigma_t = sqrt([0.01,0.1,0.5,0.9]);
L_s = length(sigma_t);
L_R = 10000;
rou = 0.95;  %%Э����������ɵĳ�������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
R = fun_rho(rou,N,1);
L=round(n*N); 
h = waitbar(0,'Please wait...');
alpha = zeros(L_s,L_R);
alpha_T = zeros(L_s,L_R);
alpha_S = zeros(L_s,L_R);
alpha_P = zeros(L_s,L_R);
alpha_ML = zeros(L_s,L_R);
for i_s = 1:L_s%%sigma
    waitbar(i_s/L_s,h,sprintf([num2str(i_s/L_s*100),'%%']));
    parfor i =1:L_R
        R_KA = zeros(size(R));
        t = normrnd(1,sigma_t(i_s),N,1);%%0~0.5%%ʧ������
        R_KA = R.*(t*t');
        Train = fun_TrainData(str_train,N,L,R,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
        x0 = fun_TrainData(str_train,N,1,R,lambda,mu,opt_train); 
        R_SCM = (fun_SCMN(Train));
        [R_CC,alpha(i_s,i)]=fun_CC(Train,R_SCM,R_KA);
        [R_ECCT,alpha_T(i_s,i)] = fun_PowerCC(Train,R_KA,1,4);
        [R_ECCS,alpha_S(i_s,i)] = fun_PowerCC(Train,R_KA,1,8);
        [R_ECCP,alpha_P(i_s,i)] = fun_PowerCC(Train,R_KA,1,7);
        [R_ML,alpha_ML(i_s,i)]=fun_MLalpha(Train,R_SCM,R_KA,x0);
    end
end
close(h)
for i = 1:L_s
    figure(i)
    hold on
    subplot(5,1,1)
    plot(alpha(i,:),'b','LineWidth',1)%k-.
    ylabel('\alpha','FontSize',10)
    axis([1,L_R,0,1])
    h_leg=legend('CC');
    axis([1,L_R,0,1])
    set(h_leg,'Location','NorthEast')
    grid on
    box on
    
    subplot(5,1,2)
    plot(alpha_ML(i,:),'c','LineWidth',1)%k:
    ylabel('\alpha','FontSize',10)
    axis([1,L_R,0,1])
    h_leg=legend('ML');
    set(h_leg,'Location','NorthEast')
    grid on
    box on
    
    subplot(5,1,3)
    plot(alpha_T(i,:),'r','LineWidth',1)
    ylabel('\alpha','FontSize',10)
    axis([1,L_R,0,1])
    h_leg=legend('T-KA');
    set(h_leg,'Location','NorthEast')
    grid on
    box on
    
    subplot(5,1,4)
    plot(alpha_S(i,:),'y','LineWidth',1)
    ylabel('\alpha','FontSize',10)
    axis([1,L_R,0,1])
    h_leg=legend('S-KA');
    set(h_leg,'Location','NorthEast')
    grid on
    box on
    
    subplot(5,1,5)
    plot(alpha_P(i,:),'g','LineWidth',1)
    ylabel('\alpha','FontSize',10)
    h_leg=legend('KA-P');
    set(h_leg,'Location','NorthEast')
    grid on
    box on
    xlabel('Range Index','FontSize',10)
end
% str = [str_train,'_','Alpha_each','_',num2str(n),'N','.mat'];
% save (str); 


