clc
clear
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
load p_fc_SNR_15_1N_s0.1.mat
figure(1);
hold on
plot(fc,Pd_NMF_mc,'r','linewidth',2)
plot(fc,Pd_SCM_mc,'g','linewidth',2)
plot(fc,Pd_CC_mc,'b','linewidth',2)
plot(fc,Pd_ML_mc,'c','linewidth',2)
plot(fc,Pd_ECCT_mc,'k','linewidth',2)
plot(fc,Pd_ECCS_mc,'K-*','linewidth',2)
plot(fc,Pd_ECCP_mc,'k-o','linewidth',2)
h_leg = legend('NMF','ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCT','ANMF with ECCS','ANMF with ECCP');
xlabel('f_{c}','FontSize',20)
ylabel('Pd','FontSize',20)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[100 0 1200 1000])
set(h_leg,'Location','NorthWest')
grid on
box on