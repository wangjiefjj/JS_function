clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 10;
FontSize = 20;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%-5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
% load ROC_IPIX4Second.mat
% plot(PFA(3:end),Pd_ECC_Mlti_mc((3:end),1),'r-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_PCC_Mlti_mc((3:end),1),'g-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_LogCC_Mlti_mc((3:end),1),'b-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_E_Mlti_mc((3:end),1),'r-.','linewidth',linewide1,'markersize',mkft)
% % load ROC_LogMP4Second_p.mat
% plot(PFA(3:end),Pd_P_Mlti_mc((3:end),1),'g-.','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_LogM_Mlti_mc((3:end),1),'b-.','linewidth',linewide1,'markersize',mkft)
% % load ROC_4Second_s0.1_p.mat
% plot(PFA(3:end),Pd_CC_Mlti_mc((3:end),1),'k-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_SFP_Mlti_mc((3:end),1),'c-*','linewidth',linewide1,'markersize',mkft)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(PFA(3:end),Pd_ECC_Mlti_mc((3:end),2),'r-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_PCC_Mlti_mc((3:end),2),'g-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_LogCC_Mlti_mc((3:end),2),'b-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_E_Mlti_mc((3:end),2),'r-.','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_P_Mlti_mc((3:end),2),'g-.','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_LogM_Mlti_mc((3:end),2),'b-.','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_CC_Mlti_mc((3:end),2),'k-*','linewidth',linewide1,'markersize',mkft)
% plot(PFA(3:end),Pd_SFP_Mlti_mc((3:end),2),'c-*','linewidth',linewide1,'markersize',mkft)
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%10dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot(PFA,Pd_ECC_Mlti_mc(:,3),'r-*','linewidth',linewide1,'markersize',mkft)
% % plot(PFA,Pd_PCC_Mlti_mc(:,3),'g-*','linewidth',linewide1,'markersize',mkft)
% % plot(PFA,Pd_LogCC_Mlti_mc(:,3),'b-*','linewidth',linewide1,'markersize',mkft)
% % plot(PFA,Pd_E_Mlti_mc(:,3),'r-.','linewidth',linewide1,'markersize',mkft)
% % plot(PFA,Pd_P_Mlti_mc(:,3),'g-.','linewidth',linewide1,'markersize',mkft)
% % plot(PFA,Pd_LogM_Mlti_mc(:,3),'b-.','linewidth',linewide1,'markersize',mkft)
% % plot(PFA,Pd_CC_Mlti_mc(:,3),'k-*','linewidth',linewide1,'markersize',mkft)
% % plot(PFA,Pd_SFP_Mlti_mc(:,3),'c-*','linewidth',linewide1,'markersize',mkft)
% % xlabel('P_{fa}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% % ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[700 0 1200 1000])
% % axis([-5,20,0,1])
% grid on
% grid minor
% box on
% xlabel('P_{fa}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% axis([min(PFA(3:end)),max(PFA(3:end)),0,1])
% h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
%     'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
%     'ANMF with SFP','NMF');
% set(h_leg,'Location','SouthEast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load ROC_IPIX_PFA_4Second.mat
plot(log10(PFA),Pd_ECC_Mlti_mc(:,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_PCC_Mlti_mc(:,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogCC_Mlti_mc(:,1),'b-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_E_Mlti_mc(:,1),'r-.','linewidth',linewide1,'markersize',mkft)
% load ROC_LogMP4Second_p.mat
plot(log10(PFA),Pd_P_Mlti_mc(:,1),'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogM_Mlti_mc(:,1),'b-.','linewidth',linewide1,'markersize',mkft)
% load ROC_4Second_s0.1_p.mat
plot(log10(PFA),Pd_CC_Mlti_mc(:,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_SFP_Mlti_mc(:,1),'c-*','linewidth',linewide1,'markersize',mkft)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(log10(PFA),Pd_ECC_Mlti_mc(:,2),'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_PCC_Mlti_mc(:,2),'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogCC_Mlti_mc(:,2),'b-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_E_Mlti_mc(:,2),'r-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_P_Mlti_mc(:,2),'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_LogM_Mlti_mc(:,2),'b-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_CC_Mlti_mc(:,2),'k-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Pd_SFP_Mlti_mc(:,2),'c-*','linewidth',linewide1,'markersize',mkft)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%10dB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(log10(PFA),Pd_ECC_Mlti_mc(:,3),'r-*','linewidth',linewide1,'markersize',mkft)
% plot(log10(PFA),Pd_PCC_Mlti_mc(:,3),'g-*','linewidth',linewide1,'markersize',mkft)
% plot(log10(PFA),Pd_LogCC_Mlti_mc(:,3),'b-*','linewidth',linewide1,'markersize',mkft)
% plot(log10(PFA),Pd_E_Mlti_mc(:,3),'r-.','linewidth',linewide1,'markersize',mkft)
% plot(log10(PFA),Pd_P_Mlti_mc(:,3),'g-.','linewidth',linewide1,'markersize',mkft)
% plot(log10(PFA),Pd_LogM_Mlti_mc(:,3),'b-.','linewidth',linewide1,'markersize',mkft)
% plot(log10(PFA),Pd_CC_Mlti_mc(:,3),'k-*','linewidth',linewide1,'markersize',mkft)
% plot(log10(PFA),Pd_SFP_Mlti_mc(:,3),'c-*','linewidth',linewide1,'markersize',mkft)
% xlabel('P_{fa}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
xlabel('Log_{PFA}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(log10(PFA)),max(log10(PFA)),0,1])
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthEast')
