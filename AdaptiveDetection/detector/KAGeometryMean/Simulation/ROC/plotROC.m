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
load ROC_4Second_s0.1_p.mat
plot(PFA,Pd_ECC_Mlti_mc(:,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_PCC_Mlti_mc(:,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_LogCC_Mlti_mc(:,1),'b-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_E_Mlti_mc(:,1),'r-.','linewidth',linewide1,'markersize',mkft)
% load ROC_LogMP4Second_p.mat
plot(PFA,Pd_P_Mlti_mc(:,1),'g-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_LogM_Mlti_mc(:,1),'b-.','linewidth',linewide1,'markersize',mkft)
% load ROC_4Second_s0.1_p.mat
plot(PFA,Pd_CC_Mlti_mc(:,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_SFP_Mlti_mc(:,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_R_Mlti_mc(:,1),'m-.*','linewidth',linewide1,'markersize',mkft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.1,5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
% load ROC_4Second_s0.1_p.mat
% Pd_P_Mlti_mc(1:2,2) = Pd_P_Mlti_mc(1:2,2)+0.2;
plot(PFA,Pd_ECC_Mlti_mc(:,2),'r-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_PCC_Mlti_mc(:,2),'g-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_LogCC_Mlti_mc(:,2),'b-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_E_Mlti_mc(:,2),'r-.','linewidth',linewide1,'markersize',mkft)
% load ROC_LogMP4Second_p.mat
plot(PFA,Pd_P_Mlti_mc(:,2),'g-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_LogM_Mlti_mc(:,2),'b-.','linewidth',linewide1,'markersize',mkft)
% load ROC_4Second_s0.1_p.mat
plot(PFA,Pd_CC_Mlti_mc(:,2),'k-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_SFP_Mlti_mc(:,2),'c-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_R_Mlti_mc(:,2),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('P_{fa}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9  -5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load ROC_4Second_s0.9_p.mat
plot(PFA,Pd_ECC_Mlti_mc(:,1),'r-*','linewidth',linewide1,'markersize',mkft)
load ROC_PCC4Second_s0.9_p.mat
plot(PFA,Pd_PCC_Mlti_mc(:,1),'g-*','linewidth',linewide1,'markersize',mkft)
load ROC_LogCC4Second_s0.9_p.mat
plot(PFA,Pd_LogCC_Mlti_mc(:,1),'b-*','linewidth',linewide1,'markersize',mkft)
load ROC_4Second_s0.9_p.mat
plot(PFA,Pd_E_Mlti_mc(:,1),'r-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_P_Mlti_mc(:,1),'g-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_LogM_Mlti_mc(:,1),'b-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_CC_Mlti_mc(:,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_SFP_Mlti_mc(:,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_R_Mlti_mc(:,1),'m-.*','linewidth',linewide1,'markersize',mkft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9,5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
load ROC_4Second_s0.9_p.mat
plot(PFA,Pd_ECC_Mlti_mc(:,2),'r-*','linewidth',linewide1,'markersize',mkft)
load ROC_PCC4Second_s0.9_p.mat
plot(PFA,Pd_PCC_Mlti_mc(:,2),'g-*','linewidth',linewide1,'markersize',mkft)
load ROC_LogCC4Second_s0.9_p.mat
plot(PFA,Pd_LogCC_Mlti_mc(:,2),'b-*','linewidth',linewide1,'markersize',mkft)
load ROC_4Second_s0.9_p.mat
plot(PFA,Pd_E_Mlti_mc(:,2),'r-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_P_Mlti_mc(:,2),'g-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_LogM_Mlti_mc(:,2),'b-.','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_CC_Mlti_mc(:,2),'k-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_SFP_Mlti_mc(:,2),'c-*','linewidth',linewide1,'markersize',mkft)
plot(PFA,Pd_R_Mlti_mc(:,2),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('P_{fa}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')