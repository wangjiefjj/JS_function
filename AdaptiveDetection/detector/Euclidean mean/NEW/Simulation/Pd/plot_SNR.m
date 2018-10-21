clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load p_1N_s0.01.mat
% figure(1);
% hold on
% plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
%     'ANMF with P-KA','ANMF with S-KA','ANMF with T-KA');
% xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[300 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load p_1N_s0.1.mat
% figure(2);
% hold on
% plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
%     'ANMF with P-KA','ANMF with S-KA','ANMF with T-KA');
% xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[300 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load p_1N_s0.5.mat
% figure(3);
% hold on
% plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
%     'ANMF with P-KA','ANMF with S-KA','ANMF with T-KA');
% xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[300 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load p_1N_s0.9.mat
% figure(4);
% hold on
% plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
%     'ANMF with P-KA','ANMF with S-KA','ANMF with T-KA');
% xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[300 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.01.mat
figure(1);
hold on
plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCP','ANMF with ECCS','ANMF with ECCT');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.1.mat
figure(2);
hold on
plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with P-KA','ANMF with S-KA','ANMF with T-KA');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.5.mat
figure(3);
hold on
plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCP','ANMF with ECCS','ANMF with ECCT');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_2N_s0.9.mat
figure(4);
hold on
plot(SNRout,Pd_NMF_mc,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCP','ANMF with ECCS','ANMF with ECCT');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on