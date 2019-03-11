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
load TSV_4Second0SNR_s0.1_p.mat
plot(TSV_ECC,'r-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_PCC,'g-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_LogCC,'b-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_E,'r-.','linewidth',linewide1,'markersize',mkft)
plot(TSV_P,'g-.','linewidth',linewide1,'markersize',mkft)
plot(TSV_LogM,'b-.','linewidth',linewide1,'markersize',mkft)
plot(TSV_CC,'k-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_SFP,'c-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_NSCM,'m-.>','linewidth',linewide1,'markersize',mkft)
plot(TSV,'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('Rang Cell','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Statistic value','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E, K=0.5N','ANMF with KA-PE, K=0.5N','ANMF with KA-LogE, K=0.5N',...
    'ANMF with E, K=0.5N','ANMF with P, K=0.5N','ANMF with LogE, K=0.5N','ANMF with CC, K=0.5N',...
    'ANMF with SFP, K=0.5N','ANMF with NSCM, K=0.5N','NMF');
set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load TSV_4Second0SNR_s0.9_p.mat

plot(TSV_ECC,'r-*','linewidth',linewide1,'markersize',mkft)%/max(TSV_ECC)
plot(TSV_PCC,'g-*','linewidth',linewide1,'markersize',mkft)%/max(TSV_PCC)
plot(TSV_LogCC,'b-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_E,'r-.','linewidth',linewide1,'markersize',mkft)
plot(TSV_P,'g-.','linewidth',linewide1,'markersize',mkft)
plot(TSV_LogM,'b-.','linewidth',linewide1,'markersize',mkft)
plot(TSV_CC,'k-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_SFP,'c-*','linewidth',linewide1,'markersize',mkft)
plot(TSV_NSCM,'m-.>','linewidth',linewide1,'markersize',mkft)
plot(TSV,'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('Rang Cell','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Statistic value','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E, K=0.5N','ANMF with KA-PE, K=0.5N','ANMF with KA-LogE, K=0.5N',...
    'ANMF with E, K=0.5N','ANMF with P, K=0.5N','ANMF with LogE, K=0.5N','ANMF with CC, K=0.5N',...
    'ANMF with SFP, K=0.5N','ANMF with NSCM, K=0.5N','NMF');
set(h_leg,'Location','NorthEast')