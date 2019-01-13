clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize2=35;
fw2 = 'normal'; %%是否加粗斜体之类
fn2='Times New Roman';
linewide2=3;
mkft2 = 10;
FontSize2 = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load PdFd_2_4Second19980223_170435_IPIX.mat
plot(fd,Pd_ECC_mc,'r-*','linewidth',linewide2,'markersize',mkft2)
plot(fd,Pd_PCC_mc,'g-*','linewidth',linewide2,'markersize',mkft2)
plot(fd,Pd_LogCC_mc,'b-*','linewidth',linewide2,'markersize',mkft2)
plot(fd,Pd_E_mc,'r-.','linewidth',linewide2,'markersize',mkft2)
plot(fd,Pd_P_mc,'g-.','linewidth',linewide2,'markersize',mkft2)
plot(fd,Pd_LogM_mc,'b-.','linewidth',linewide2,'markersize',mkft2)
plot(fd,Pd_CC_mc,'k-*','linewidth',linewide2,'markersize',mkft2)
plot(fd,Pd_SFP_mc,'c-*','linewidth',linewide2,'markersize',mkft2)

% xlabel('Normalized Doppler of Target','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
xlabel('目标归一化多普勒','FontSize',labeltsize2,'FontWeight',fw2)
ylabel('PD','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
set(gca,'FontSize',labeltsize2)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP');
set(h_leg,'Location','SouthWest')
