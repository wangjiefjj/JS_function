clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%Homogenous%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%辅助数据
load Homogenous_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Homogenous_Accuracy_4.mat
plot(L,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['SecondaryOnly'];
%%
%%主辅数据
load Homogenous_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC2,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc2,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC2,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC2,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Homogenous_Accuracy_4.mat
plot(L,Accuracy_GIC2,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PrimarySecondary.eps'];
%%
%%主数据
load Homogenous_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC3,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc3,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC3,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC3,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Homogenous_Accuracy_4.mat
plot(L,Accuracy_GIC3,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PrimaryOnly.eps'];
%%%%%%%%%%%%%%%%%%%%%Partial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%副数据
load Partial_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_Accuracy_4.mat
plot(L,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['SecondaryOnly.eps'];
%%
%%主副数据
load Partial_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC2,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc2,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC2,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC2,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_Accuracy_4.mat
plot(L,Accuracy_GIC2,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PrimarySecondary.eps'];
%%
%%主数据
load Partial_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC3,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc3,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC3,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC3,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_Accuracy_4.mat
plot(L,Accuracy_GIC3,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PrimaryOnly.eps'];
%%%%%%%%%%%%%%%%%%%%%SIRP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%辅助数据
load SIRP_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_Accuracy_4.mat
plot(L,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['SecondaryOnly.eps'];
%%
%%主副数据
load SIRP_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC2,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc2,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC2,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC2,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_Accuracy_4.mat
plot(L,Accuracy_GIC2,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PrimarySecondary.eps'];
%%
%%主数据
load SIRP_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC3,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc3,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC3,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC3,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_Accuracy_4.mat
plot(L,Accuracy_GIC3,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PrimaryOnly.eps'];