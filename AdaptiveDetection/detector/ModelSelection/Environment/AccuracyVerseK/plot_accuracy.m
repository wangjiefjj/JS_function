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
% %%辅助数据
% load Homogenous_Accuracy_2.mat
% figure();
% hold on
% plot(L,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load Homogenous_Accuracy_4.mat
% plot(L,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
% xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([9,32,0,100])
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H1_SecondaryOnly','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主辅数据
load Hom_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Hom_Accuracy_4.mat
plot(L,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{1,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([9,32,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主数据
load Hom_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Hom_Accuracy_4.mat
plot(L,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{1,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([9,32,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%Partial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%副数据
% load Partial_Accuracy_2.mat
% figure();
% hold on
% plot(L,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load Partial_Accuracy_4.mat
% plot(L,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([9,32,0,100])
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H2_SecondaryOnly.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主副数据
load Partial_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_Accuracy_4.mat
plot(L,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([9,32,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主数据
load Partial_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_Accuracy_4.mat
plot(L,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([9,32,0,100])
axis([9,32,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%SIRP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%辅助数据
% load SIR\it{P}_Accuracy_2.mat
% figure();
% hold on
% plot(L,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(L,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load SIR\it{P}_Accuracy_4.mat
% plot(L,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Recognition Percentage','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([9,32,0,100])
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H3_SecondaryOnly.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主副数据
load SIRP_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_Accuracy_4.mat
plot(L,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{3,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([9,32,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主数据
load SIRP_Accuracy_2.mat
figure();
hold on
plot(L,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(L,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_Accuracy_4.mat
plot(L,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','ABIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('K','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{3,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([9,32,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径