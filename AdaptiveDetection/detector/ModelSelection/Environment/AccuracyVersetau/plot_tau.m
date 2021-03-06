clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%H1vsH2%%%%%%%%%%%%%%%%%%%%%%%
%% 主副数据
load Hom_tau_L20_rho2.mat
figure();
hold on
plot(tau,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Hom_tau_L20_rho4.mat
plot(tau,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{2,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,10,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_tau_PrimarySecondary','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%% 主数据
% load Partial_tau_L20_rho2.mat
% figure();
% hold on
% plot(tau,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load Partial_tau_L20_rho4.mat
% plot(tau,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('\it{P}_{2,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([1,10,0,100])
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H1_tau_PrimaryOnly','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%H2vsH2%%%%%%%%%%%%%%%%%%%%%%%
%% 主副数据
load Partial_tau_L20_rho2.mat
figure();
hold on
plot(tau,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load Partial_tau_L20_rho4.mat
plot(tau,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,10,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_tau_PrimarySecondary','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%% 主数据
% load Partial_tau_L20_rho2.mat
% figure();
% hold on
% plot(tau,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load Partial_tau_L20_rho4.mat
% plot(tau,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([1,10,0,100])
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H2_tau_PrimaryOnly','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%H3vsH2%%%%%%%%%%%%%%%%%%%%%%%
%% 主副数据
load SIRP_tau_L20_rho2.mat
figure();
hold on
plot(tau,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(tau,Accuracy_GIC2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
load SIRP_tau_L20_rho4.mat
plot(tau,Accuracy_GIC2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\it{P}_{2,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([1,10,0,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_tau_PrimarySecondary','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%% 主数据
% load SIRP_tau_L20_rho2.mat
% figure();
% hold on
% plot(tau,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(tau,Accuracy_GIC3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load SIRP_tau_L20_rho4.mat
% plot(tau,Accuracy_GIC3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','wBIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('\tau','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('\it{P}_{2,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([1,10,0,100])
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H3_tau_PrimaryOnly','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径