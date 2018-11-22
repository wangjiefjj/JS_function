clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%
%%%%%%%%%%%%%%%H1VsH3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 主副数据
load Hom_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Hom_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\lambda','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_lambda_PrimarySecondary','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%% 主数据
load Hom_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Hom_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\lambda','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_lambda_PrimaryOnly','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%%%%%%%%%%%%%%H2VsH3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 主副数据
load Partial_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\lambda','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_lambda_PrimarySecondary','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%% 主数据
load Partial_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\lambda','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_lambda_PrimaryOnly','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%%%%%%%%%%%%%%H3VsH3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 主副数据
load SIRP_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\lambda','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_lambda_PrimarySecondary','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%% 主数据
load SIRP_lambda_L20_rho2.mat
figure();
hold on
plot(lambda,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(lambda,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_lambda_L20_rho4.mat
plot(lambda,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('\lambda','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,3,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_lambda_PrimaryOnly','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径