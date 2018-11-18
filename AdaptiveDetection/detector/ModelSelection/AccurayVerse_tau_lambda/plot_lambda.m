clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%�Ƿ�Ӵ�б��֮��
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%
%% ��������
load SIRP_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC2,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc2,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC2,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC2,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC2,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,1])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['SIRP_lambda_PrimarySecondary','.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%% ������
load SIRP_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC3,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc3,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC3,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC3,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC3,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,1])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['SIRP_lambda_PrimaryOnly','.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��