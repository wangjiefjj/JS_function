clc
clear 
close all
labeltsize = 35;
fw = 'normal'; %%�Ƿ�Ӵ�б��֮��
fn = 'Times New Roman';
linewide1 = 3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%Homogenous%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%��������
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
% print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%��������
load Hom_fd_L20.mat
figure();
hold on
plot(theta_sig,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC2_2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC4_2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('Ŀ���һ��Ƶ��','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{1,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(theta_sig),max(theta_sig),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_fd_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%������
load Hom_fd_L20.mat
figure();
hold on
plot(theta_sig,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC2_3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC4_3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('Ŀ���һ��Ƶ��','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{1,1}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(theta_sig),max(theta_sig),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H1_fd_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%%%%%%%%%%%%%%%%%%%%Partial%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%������
% load Partial_SNR_L20_theta_sig2.mat
% figure();
% hold on
% plot(SCNRout,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SCNRout,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load Partial_SNR_L20_theta_sig4.mat
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
% print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%��������
load Partial_fd_L20.mat
figure();
hold on
plot(theta_sig,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC2_2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC4_2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('Ŀ���һ��Ƶ��','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(theta_sig),max(theta_sig),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_fd_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%������
load Partial_fd_L20.mat
figure();
hold on
plot(theta_sig,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC2_3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC4_3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('Ŀ���һ��Ƶ��','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{2,2}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(theta_sig),max(theta_sig),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H2_fd_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%%%%%%%%%%%%%%%%%%%%SIRP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%��������
% load SIRP_fd_L20_theta_sig2.mat
% figure();
% hold on
% plot(theta_sig,Accuracy_AIC1,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(theta_sig,Accuracy_AICc1,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(theta_sig,Accuracy_ABIC1,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(theta_sig,Accuracy_GIC1,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
% load SIRP_fd_L20_theta_sig4.mat
% plot(theta_sig,Accuracy_GIC1,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
% xlabel('Ŀ���һ��Ƶ��','FontSize',labeltsize,'FontWeight',fw)
% ylabel('Accuracy','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[500 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% str=['H3_SCNR_SecondaryOnly.eps'];
% print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%��������
load SIRP_fd_L20.mat
figure();
hold on
plot(theta_sig,Accuracy_AIC2*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_AICc2*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_ABIC2*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC2_2*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC4_2*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('Ŀ���һ��Ƶ��','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{3,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(theta_sig),max(theta_sig),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_fd_PrimarySecondary.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%������
load SIRP_fd_L20.mat
figure();
hold on
plot(theta_sig,Accuracy_AIC3*100,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_AICc3*100,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_ABIC3*100,'b-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC2_3*100,'k-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(theta_sig,Accuracy_GIC4_3*100,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('AIC','AICc','BIC','GIC(\eta=2)','GIC(\eta=4)');
xlabel('Ŀ���һ��Ƶ��','FontSize',labeltsize,'FontWeight',fw)
ylabel('\it{P}_{3,3}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(theta_sig),max(theta_sig),60,100])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['H3_fd_PrimaryOnly.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��