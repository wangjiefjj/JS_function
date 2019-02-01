%%ѡ�����ݵ����ģ��
%%K����,������SCNR�ı仯
clc
clear
close all
Class=3; %%
MC = 10000;
rho = 0.9;  %%Э����������ɵĳ�������
fc = 0.1;
%%%%�����������
Na = 2;     % ��Ԫ��
Np = 4;     % ������
N = Na*Np;
lambda = 1;%%%ԽС�Ǹ�˹Խ����
mu = 1;
% n =1.1:0.3:8; %����������
% L=round(n*N);SNRout=10;
SCNRout=0;
SCNRnum=10.^(SCNRout/10);
CNRout=30;
CNRnum=10.^(CNRout/10); 
L = 20;
theta_sig = -0.5:0.05:0.5;
nn = 0:N-1;
if Class==1%%����
    str_train = 'g';
     %%�Ӳ�Э����
    Rc1 = fun_rho(rho,N,1,fc);
    %     sigmaf = 0.03; %%�Ӳ���չ��
    %     rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
    %     Rc1 = CNRnum * toeplitz(rc);
    Rc1 = Rc1+ 1/CNRnum * eye(N) ;%
    Rc2 = Rc1;
    opt_train = 0;
    str=['Hom_fd_L',num2str(L),'.mat'];
    elseif Class == 2%%���־���
        str_train = 'g';
        opt_train = 2;
     %%�Ӳ�Э����
        Rc2 = fun_rho(rho,N,1,fc);
    %     sigmaf = 0.03; %%�Ӳ���չ��
    %     rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
    %     Rc1 = CNRnum * toeplitz(rc);
        Rc2 = Rc2+ 1/CNRnum * eye(N) ;%+ eye(N)
        Rc1 = 10*Rc2;
        str=['Partial_fd_L',num2str(L),'.mat'];
    elseif Class == 3%%SIRP
        str_train = 'p';
     %%%IG��ѡ�1Ϊÿ�����뵥ԪIG������ͬ
     %2ÿ����Ԫ����ֵһ����Ϊ���־��Ȼ���
        opt_train = 1;     
     %%�Ӳ�Э����
        Rc1 = fun_rho(rho,N,1,fc);
    %     sigmaf = 0.03; %%�Ӳ���չ��
    %     rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
    %     Rc1 = CNRnum * toeplitz(rc);
        Rc1 = Rc1+ 1/CNRnum * eye(N) ;%+ eye(N)
        Rc2 = Rc1;
        str=['SIRP_fd_L',num2str(L),'.mat'];
end
h = waitbar(1,'Please wait...');
tic
for i_fd = 1:length(theta_sig)
    waitbar(i_fd/length(theta_sig),h,sprintf([num2str(i_fd/length(theta_sig)*100),'%%']));
    p = exp(-1i*2*pi*nn*theta_sig(i_fd)).'/sqrt(N); %%%%%% ϵͳ����ʸ��
    iRc1 = inv(Rc1);
    iRc2 = inv(Rc2);
    a = sqrt(SCNRnum/abs(p'/Rc2*p));
    count_AIC1 = 0;
    count_ABIC1 = 0;
    count_GIC2_1 = 0;
    count_GIC4_1 = 0;
    count_AICc1 = 0;

    count_AIC2 = 0;
    count_ABIC2 = 0;
    count_GIC2_2 = 0;
    count_GIC4_2 = 0;
    count_AICc2 = 0;

    count_AIC3 = 0;
    count_ABIC3 = 0;
    count_GIC2_3 = 0;
    count_GIC4_3 = 0;
    count_AICc3 = 0;
     %%��������
     %ֻ�ø�������ʱ
    H1_num1 = N^2;
    H2_num1 = N^2+1;
    H3_num1 = N^2+L;
     %��������ʱ
    H1_num2 = N^2+1;
    H2_num2 = N^2+3;
    H3_num2 = N^2+L+2;
     %������ʱ
    H1_num3 = N^2+1;
    H2_num3 = N^2+2;
    H3_num3 = N^2+2;
   parfor i = 1:MC
      warning off
      %%��������
      [Train,tauk] = fun_TrainData(str_train,N,L,Rc1,lambda,mu,opt_train);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
      [x0,tau0] = fun_TrainData(str_train,N,1,Rc2,lambda,mu,opt_train); % �����źŽ������Ӳ�������
      x0 =x0+a*p;
      %%s�������� 
      %ֻ�ø�������ʱ
      s_H1_1 = -abs(fun_s_H1(Train,p,1));
      s_H2_1 = -abs(fun_s_H2(Train,p,1));
      [s_H3_1] = -abs(fun_s_H3(Train,p,1));
      %��������ʱ
      s_H1_2 = -abs(fun_s_H1([Train,x0],p,2));
      s_H2_2 = -abs(fun_s_H2([Train,x0],p,2));
      s_H3_2 = -abs(fun_s_H3([Train,x0],p,2));
      %ֻ��������ʱ
      s_H1_3 = -abs(fun_s_H1([Train,x0],p,3));
      s_H2_3 = -abs(fun_s_H2([Train,x0],p,3));
      s_H3_3 = -abs(fun_s_H3([Train,x0],p,3));

      %%ģ��ѡ��׼��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %ֻ�ø�������ʱ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%AIC%%%%%%%%%
      H1_AIC_1 = fun_AIC(s_H1_1,H1_num1);
      H2_AIC_1 = fun_AIC(s_H2_1,H2_num1);
      H3_AIC_1 = fun_AIC(s_H3_1,H3_num1);
      Class_AIC1 =[H1_AIC_1,H2_AIC_1,H3_AIC_1];
      [~,Class_AIC_num1] = min(Class_AIC1);
      if Class_AIC_num1 == Class  
          count_AIC1 = count_AIC1+1;  
      end
      %%GIC2%%%%%%%%%%%%%%
      H1_GIC2_1 = fun_GIC(s_H1_1,H1_num1,2);
      H2_GIC2_1 = fun_GIC(s_H2_1,H2_num1,2);
      H3_GIC2_1 = fun_GIC(s_H3_1,H3_num1,2);
      Class_GIC2_1 =[H1_GIC2_1,H2_GIC2_1,H3_GIC2_1];
      [~,Class_GIC2_num1] = min(Class_GIC2_1);
      if Class_GIC2_num1 == Class  
          count_GIC2_1 = count_GIC2_1+1;  
      end
      %%GIC4%%%%%%%%%%%%%%
      H1_GIC4_1 = fun_GIC(s_H1_1,H1_num1,4);
      H2_GIC4_1 = fun_GIC(s_H2_1,H2_num1,4);
      H3_GIC4_1 = fun_GIC(s_H3_1,H3_num1,4);
      Class_GIC4_1 =[H1_GIC4_1,H2_GIC4_1,H3_GIC4_1];
      [~,Class_GIC4_num1] = min(Class_GIC4_1);
      if Class_GIC4_num1 == Class  
          count_GIC4_1 = count_GIC4_1+1;  
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
      %��������ʱ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%AIC%%%%%%%%%%%%%%%%
      AIC_H1_2 = fun_AIC(s_H1_2,H1_num2);
      AIC_H2_2 = fun_AIC(s_H2_2,H2_num2);
      AIC_H3_2 = fun_AIC(s_H3_2,H3_num2);
      Class_AIC2 =[AIC_H1_2,AIC_H2_2,AIC_H3_2];
      [~,Class_AIC_num2] = min(Class_AIC2);
      if Class_AIC_num2 == Class  
          count_AIC2 = count_AIC2+1;  
      end
     %%GIC2%%%%%%%%%%%%%%
      H1_GIC2_2 = fun_GIC(s_H1_2,H1_num2,2);
      H2_GIC2_2 = fun_GIC(s_H2_2,H2_num2,2);
      H3_GIC2_2 = fun_GIC(s_H3_2,H3_num2,2);
      Class_GIC2_2 =[H1_GIC2_2,H2_GIC2_2,H3_GIC2_2];
      [~,Class_GIC2_num2] = min(Class_GIC2_2);
      if Class_GIC2_num2 == Class  
          count_GIC2_2 = count_GIC2_2+1;  
      end
      %%GIC4%%%%%%%%%%%%%%
      H1_GIC4_2 = fun_GIC(s_H1_2,H1_num2,4);
      H2_GIC4_2 = fun_GIC(s_H2_2,H2_num2,4);
      H3_GIC4_2 = fun_GIC(s_H3_2,H3_num2,4);
      Class_GIC4_2 =[H1_GIC4_2,H2_GIC4_2,H3_GIC4_2];
      [~,Class_GIC4_num2] = min(Class_GIC4_2);
      if Class_GIC4_num2 == Class  
          count_GIC4_2 = count_GIC4_2+1;  
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
      %������ʱ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%AIC%%%%%%%%%%%%%%%%
      AIC_H1_3 = fun_AIC(s_H1_3,H1_num3);
      AIC_H2_3 = fun_AIC(s_H2_3,H2_num3);
      AIC_H3_3 = fun_AIC(s_H3_3,H3_num3);
      Class_AIC3 =[AIC_H1_3,AIC_H2_3,AIC_H3_3];
      [~,Class_AIC_num3] = min(Class_AIC3);
      if Class_AIC_num3 == Class  
          count_AIC3 = count_AIC3+1;  
      end
      %%GIC2%%%%%%%%%%%%%%
      H1_GIC2_3 = fun_GIC(s_H1_3,H1_num3,2);
      H2_GIC2_3 = fun_GIC(s_H2_3,H2_num3,2);
      H3_GIC2_3 = fun_GIC(s_H3_3,H3_num3,2);
      Class_GIC2_3 =[H1_GIC2_3,H2_GIC2_3,H3_GIC2_3];
      [~,Class_GIC2_num3] = min(Class_GIC2_3);
      if Class_GIC2_num3 == Class  
          count_GIC2_3 = count_GIC2_3+1;  
      end
      %%GIC4%%%%%%%%%%%%%%
      H1_GIC4_3 = fun_GIC(s_H1_3,H1_num3,4);
      H2_GIC4_3 = fun_GIC(s_H2_3,H2_num3,4);
      H3_GIC4_3 = fun_GIC(s_H3_3,H3_num3,4);
      Class_GIC4_3 =[H1_GIC4_3,H2_GIC4_3,H3_GIC4_3];
      [~,Class_GIC4_num3] = min(Class_GIC4_3);
      if Class_GIC4_num3 == Class  
          count_GIC4_3 = count_GIC4_3+1;  
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
     Accuracy_AIC1(i_fd) = count_AIC1/MC;
     Accuracy_GIC2_1(i_fd) = count_GIC2_1/MC;
     Accuracy_GIC4_1(i_fd) = count_GIC4_1/MC;
     Accuracy_AICc1(i_fd) = count_AICc1/MC;
     Accuracy_ABIC1(i_fd) = count_ABIC1/MC;

     Accuracy_AIC2(i_fd) = count_AIC2/MC;
     Accuracy_GIC2_2(i_fd) = count_GIC2_2/MC;
     Accuracy_GIC4_2(i_fd) = count_GIC4_2/MC;
     Accuracy_AICc2(i_fd) = count_AICc2/MC;
     Accuracy_ABIC2(i_fd) = count_ABIC2/MC;

     Accuracy_AIC3(i_fd) = count_AIC3/MC;
     Accuracy_GIC2_3(i_fd) = count_GIC2_3/MC;
     Accuracy_GIC4_3(i_fd) = count_GIC4_3/MC;
     Accuracy_AICc3(i_fd) = count_AICc3/MC;
     Accuracy_ABIC3(i_fd) = count_ABIC3/MC;
end
toc
close(h)
% figure(1)
% hold on
% plot(rho,Accuracy_AIC1,'r-s','LineWidth',2)
% plot(rho,Accuracy_GIC1,'g-o','LineWidth',2)
% plot(rho,Accuracy_AICc1,'b-*','LineWidth',2)
% plot(rho,Accuracy_ABIC1,'k-.','LineWidth',2)
% title('ֻ�ø�������')
% h_leg1 = legend('AIC','GIC','AICc','ABIC');
% xlabel('K')
% ylabel('Accuracy')
% set(gca,'FontSize',10)
% set(h_leg1,'Location','SouthEast')
% grid on
% box on
figure(2)
hold on
plot(theta_sig,Accuracy_AIC2,'r-s','LineWidth',2)
plot(theta_sig,Accuracy_GIC2_2,'g-o','LineWidth',2)
plot(theta_sig,Accuracy_GIC4_2,'g-.o','LineWidth',2)
plot(theta_sig,Accuracy_AICc2,'b-*','LineWidth',2)
plot(theta_sig,Accuracy_ABIC2,'k-.','LineWidth',2)
title('��������')
h_leg2 = legend('AIC','GIC2','GIC4','AICc','ABIC');
xlabel('K')
ylabel('Accuracy')
axis([min(theta_sig),max(theta_sig),0,1])
set(gca,'FontSize',10)
set(h_leg2,'Location','SouthEast')
grid on
box on

figure(3)
hold on
plot(theta_sig,Accuracy_AIC2,'r-s','LineWidth',2)
plot(theta_sig,Accuracy_GIC2_3,'g-o','LineWidth',2)
plot(theta_sig,Accuracy_GIC4_3,'g-.o','LineWidth',2)
plot(theta_sig,Accuracy_AICc2,'b-*','LineWidth',2)
plot(theta_sig,Accuracy_ABIC2,'k-.','LineWidth',2)
title('��������')
h_leg2 = legend('AIC','GIC2','GIC4','AICc','ABIC');
xlabel('K')
ylabel('Accuracy')
axis([min(theta_sig),max(theta_sig),0,1])
set(gca,'FontSize',10)
set(h_leg2,'Location','SouthEast')
grid on

grid on
box on
save(str)
