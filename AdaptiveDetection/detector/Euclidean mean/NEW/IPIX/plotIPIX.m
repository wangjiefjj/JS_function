clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%实测数据IPIX%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;

figure(1)
hold on
load IPIX_19980223_170435_1N.mat
plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCP','ANMF with ECCS','ANMF with ECCT');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load IPIX_19980223_170435_2N.mat
plot(SNRout,Pd_SCM_mc,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with ECCP','ANMF with ECCS','ANMF with ECCT');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on