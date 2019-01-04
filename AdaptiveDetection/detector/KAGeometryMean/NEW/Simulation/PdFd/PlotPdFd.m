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
load PdFd_4Second_p.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on

load PdFd_2_4Second_s0.1_p.mat
% Pd_E_mc(11) = Pd_E_mc(11)-0.1;
plot(fd,circshift(Pd_ECC_mc,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_PCC_mc,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_LogCC_mc,1),'b-*','linewidth',linewide1,'markersize',mkft)
load PdFd_P_E_4Second_p.mat
plot(fd,circshift(Pd_E_mc,1),'r-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_P_mc,1),'g-.','linewidth',linewide1,'markersize',mkft)
load PdFd_2_4Second_s0.1_p.mat
plot(fd,circshift(Pd_LogM_mc,1),'b-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_CC_mc,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_SFP_mc,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_R_mc,1),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('Normalized Doppler','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load PdFd_2_4Second_s0.9_p.mat
% Pd_E_mc(11) = Pd_E_mc(11)-0.1;
% t = Pd_CC_mc(12);
% Pd_CC_mc(12) = Pd_CC_mc(11);
% Pd_CC_mc(11) = t;
% Pd_CC_mc(12:end-1) = Pd_CC_mc(12:end-1) + 0.08;
plot(fd,circshift(Pd_ECC_mc,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_PCC_mc,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_LogCC_mc,1),'b-*','linewidth',linewide1,'markersize',mkft)
load PdFd_P_E_4Second_p.mat
plot(fd,circshift(Pd_E_mc,1),'r-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_P_mc,1),'g-.','linewidth',linewide1,'markersize',mkft)
load PdFd_2_4Second_s0.9_p.mat
plot(fd,circshift(Pd_LogM_mc,1),'b-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_CC_mc,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_SFP_mc,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_R_mc,1),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('Normalized Doppler','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
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
