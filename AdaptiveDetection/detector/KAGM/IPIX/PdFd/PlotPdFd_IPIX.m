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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ϸ�˹0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load PdFd_PFA_4Second19980223_170435_IPIX.mat
plot(fd,Pd_ECC_mc,'r-*','linewidth',linewide1,'markersize',mkft)
plot(fd,Pd_PCC_mc,'g-*','linewidth',linewide1,'markersize',mkft)
plot(fd,Pd_LogCC_mc,'b-*','linewidth',linewide1,'markersize',mkft)
plot(fd,Pd_E_mc,'r-.','linewidth',linewide1,'markersize',mkft)
plot(fd,Pd_P_mc,'g-.','linewidth',linewide1,'markersize',mkft)
plot(fd,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
plot(fd,Pd_CC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
plot(fd,Pd_SFP_mc,'c-*','linewidth',linewide1,'markersize',mkft)

xlabel('Normalized Doppler of Target','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP');
set(h_leg,'Location','SouthWest')
% set(h_leg,'box','off')