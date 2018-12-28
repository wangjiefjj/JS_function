clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.01.mat
figure(1);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD1N_p_1','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.1.mat
figure(2);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD1N_p_2','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.5.mat
figure(3);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD1N_p_3','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.9.mat
figure(4);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD1N_p_4','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.01.mat
figure(5);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD2N_p_1','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.1.mat
figure(6);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD2N_p_2','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.5.mat
figure(7);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD2N_p_3','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.9.mat
figure(8);
hold on
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
grid on
box on
% str=['PD2N_p_4','.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
% close all