clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=35;
fw = 'normal'; %%�Ƿ�Ӵ�б��֮��
fn='Times New Roman';
linewide1=3;
mkft = 10;
FontSize = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ϸ�˹0.1��-5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load ROC4Second_s0.1_p.mat
plot(log10(PFA),Pd_ECC_Mlti_mc(:,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_PCC_Mlti_mc(:,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogCC_Mlti_mc(:,1),'b-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_E_Mlti_mc(:,1),'r-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_P_Mlti_mc(:,1),'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogM_Mlti_mc(:,1),'b-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_CC_Mlti_mc(:,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_SFP_Mlti_mc(:,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_R_Mlti_mc(:,1),'m-.*','linewidth',linewide1,'markersize',mkft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ϸ�˹0.1,5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(log10(PFA),Pd_ECC_Mlti_mc(:,2),'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_PCC_Mlti_mc(:,2),'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogCC_Mlti_mc(:,2),'b-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_E_Mlti_mc(:,2),'r-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_P_Mlti_mc(:,2),'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogM_Mlti_mc(:,2),'b-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_CC_Mlti_mc(:,2),'k-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_SFP_Mlti_mc(:,2),'c-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_R_Mlti_mc(:,2),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('log_{PFA}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
axis([min(log10(PFA)),max(log10(PFA)),0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ϸ�˹0.9  -5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load ROC4Second_s0.9_p.mat
plot(log10(PFA),Pd_ECC_Mlti_mc(:,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_PCC_Mlti_mc(:,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogCC_Mlti_mc(:,1),'b-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_E_Mlti_mc(:,1),'r-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_P_Mlti_mc(:,1),'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogM_Mlti_mc(:,1),'b-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_CC_Mlti_mc(:,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_SFP_Mlti_mc(:,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_R_Mlti_mc(:,1),'m-.*','linewidth',linewide1,'markersize',mkft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ϸ�˹0.9,5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
load ROC4Second_s0.9_p.mat
plot(log10(PFA),Pd_ECC_Mlti_mc(:,2),'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_PCC_Mlti_mc(:,2),'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogCC_Mlti_mc(:,2),'b-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_E_Mlti_mc(:,2),'r-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_P_Mlti_mc(:,2),'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogM_Mlti_mc(:,2),'b-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_CC_Mlti_mc(:,2),'k-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_SFP_Mlti_mc(:,2),'c-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_R_Mlti_mc(:,2),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('log_{PFA}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
axis([min(log10(PFA)),max(log10(PFA)),0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')