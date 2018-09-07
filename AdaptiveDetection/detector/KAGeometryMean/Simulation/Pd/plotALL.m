clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
FontSize = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load Pd_4Second_s0.1_p.mat
plot(SNRout,Pd_ECC_mc,'r-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_PCC_mc,'g-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_E_mc,'r-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_P_mc,'g-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_CC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_SFP_mc,'c-*','linewidth',linewide1,'markersize',mkft)
load OPT_p.mat
plot(SNRout,Pd_R_mc,'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
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
load Pd_4Second_s0.9_p.mat
plot(SNRout,Pd_ECC_mc,'r-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_PCC_mc,'g-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_E_mc,'r-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_P_mc,'g-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_CC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_SFP_mc,'c-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_NSCM_mc,'m-.','linewidth',linewide1,'markersize',mkft)
load OPT_p.mat
plot(SNRout,Pd_R_mc,'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
load Pd_10Second_s0.1_p.mat
plot(SNRout,Pd_ECC_mc,'r-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_PCC_mc,'g-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_E_mc,'r-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_P_mc,'g-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_CC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_SFP_mc,'c-*','linewidth',linewide1,'markersize',mkft)
load OPT_p.mat
plot(SNRout,Pd_R_mc,'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
hold on
load Pd_new_10Second_s0.9_p.mat

plot(SNRout(1:26),Pd_ECC_mc(2:27),'r-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout(1:26),Pd_PCC_mc(2:27),'g-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout(1:26),Pd_LogCC_mc(2:27),'b-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_E_mc,'r-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_P_mc,'g-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_CC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_SFP_mc,'c-*','linewidth',linewide1,'markersize',mkft)
load OPT_p.mat
plot(SNRout,Pd_R_mc,'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')