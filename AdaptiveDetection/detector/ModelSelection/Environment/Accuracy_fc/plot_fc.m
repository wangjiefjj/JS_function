clc
clear 
close all
labeltsize = 35;
fw = 'normal'; %%是否加粗斜体之类
fn = 'Times New Roman';
linewide1 = 3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%Homogenous%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%辅助数据
% load Hom_SNR_L20_rho2.mat
% figure();
% hold on
% plot(SCNRout,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load Hom_SNR_L20_rho4.mat
% plot(SCNRout,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
% xlabel('SCNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([-20,20,0,100])
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H1_SCNR_SecondaryOnly','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主辅数据
load Hom_fc_L20.mat
figure();
hold on
plot(fc,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC2_2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC4_2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('杂波归一化频率','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{1,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(fc),max(fc),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_fc_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主数据
load Hom_fc_L20.mat
figure();
hold on
plot(fc,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC2_3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC4_3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('杂波归一化频率','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{1,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(fc),max(fc),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_fc_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%Partial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%副数据
% load Partial_SNR_L20_fc2.mat
% figure();
% hold on
% plot(SCNRout,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load Partial_SNR_L20_fc4.mat
% plot(SCNRout,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('SCNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H2_SecondaryOnly.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主副数据
load Partial_fc_L20.mat
figure();
hold on
plot(fc,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC2_2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC4_2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('杂波归一化频率','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(fc),max(fc),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_fc_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主数据
load Partial_fc_L20.mat
figure();
hold on
plot(fc,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC2_3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC4_3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('杂波归一化频率','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(fc),max(fc),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_fc_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%SIRP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%辅助数据
% load SIRP_fc_L20_fc2.mat
% figure();
% hold on
% plot(fc,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(fc,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(fc,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(fc,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load SIRP_fc_L20_fc4.mat
% plot(fc,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('杂波归一化频率','FontSize',labeltsize,'FontWeight',fw)
% ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H3_SCNR_SecondaryOnly.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主副数据
load SIRP_fc_L20.mat
figure();
hold on
plot(fc,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC2_2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC4_2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('杂波归一化频率','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{3,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(fc),max(fc),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_fc_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%主数据
load SIRP_fc_L20.mat
figure();
hold on
plot(fc,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC2_3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(fc,Accuracy_GIC4_3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('杂波归一化频率','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{3,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(fc),max(fc),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_fc_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径