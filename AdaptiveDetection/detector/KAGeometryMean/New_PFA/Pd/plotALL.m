
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear 
close all
labeltsize2=35;
fw2 = 'normal'; %%是否加粗斜体之类
fn2='Times New Roman';
linewide2=3;
mkft2 = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load PD_8Second_s0.1_p.mat
plot(SNRout,Pd_ECC_mc,'r-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_PCC_mc,'g-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_E_mc,'r-.','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_P_mc,'g-.','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_CC_mc,'k-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_SFP_mc,'c-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_R_mc,'m-.*','linewidth',linewide2,'markersize',mkft2)

xlabel('SCR/dB','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
ylabel('PD','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
set(gca,'FontSize',labeltsize2)
set(gcf,'Position',[300 0 1200 1000])
axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load PD_4Second_s0.9_p.mat
plot(SNRout,Pd_ECC_mc,'r-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_PCC_mc,'g-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_E_mc,'r-.','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_P_mc,'g-.','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_CC_mc,'k-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_SFP_mc,'c-*','linewidth',linewide2,'markersize',mkft2)
plot(SNRout,Pd_R_mc,'m-.*','linewidth',linewide2,'markersize',mkft2)

xlabel('SCR/dB','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
ylabel('PD','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
set(gca,'FontSize',labeltsize2)
set(gcf,'Position',[300 0 1200 1000])
axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')
