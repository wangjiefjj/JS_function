clc
clear 
close all
labeltsize=10;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=1;
mkft = 5;
ha = fun_tight_subplot(2,4,[.07 .05],[.1 .01],[.03 .01]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.01.mat
set(gcf, 'PaperPositionMode', 'auto');
% subplot(2,4,1)
% subplot('Position',[0.05 0.75 0.2 0.2]);
axes(ha(1))
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
hold on
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel({'SNR/dB';'(a)'},'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gca,'LooseInset',get(gca,'TightInset'))
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,4,2,1,1);
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.1.mat
% subplot(2,4,2)
axes(ha(2))
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
hold on
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel({'SNR/dB';'(b)'},'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,2,4,1,2);
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.5.mat
% subplot(2,4,3)
axes(ha(3))
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
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,2,4,1,3);
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_1N_s0.9.mat
% subplot(2,4,4)
axes(ha(4))
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
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,2,4,1,4);
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.01.mat
% subplot(2,4,5)
axes(ha(5))
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
hold on
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,2,4,2,1);
grid on
box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.1.mat
% subplot(2,4,6)
axes(ha(6))
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
hold on
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,2,4,2,2);
grid on
box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.5.mat
% subplot(2,4,7)
axes(ha(7))
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
hold on
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,2,4,2,3);
grid on
box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.9.mat
% subplot(2,4,8)
axes(ha(8))
plot(SNRout,Pd_NMF_mc,'k','linewidth',linewide1,'MarkerSize',mkft)
hold on
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(h_leg,'Location','SouthEast')
axis([min(SNRout),max(SNRout),0,1])
% fun_RemoveSubplotWhiteArea(gca,2,4,2,4);
grid on
box on
set(gcf,'Position',[0 0 1850 1000])
% str=['PD_p','.eps'];
% print(gcf,'-deps','-r300',str)   %保存为png格式的图片到当前路径