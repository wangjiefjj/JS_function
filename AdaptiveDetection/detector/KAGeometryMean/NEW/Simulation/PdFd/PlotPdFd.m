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
% load PdFd_4Second_p.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on

load PdFd_4Second_s0.1_p.mat
plot(fd,circshift(Pd_ECC_mc,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_PCC_mc,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_LogCC_mc,1),'b-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_E_mc,1),'r-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_P_mc,1),'g-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_LogM_mc,1),'b-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_CC_mc,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_SFP_mc,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_R_mc,1),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('Normalized Doppler of Target','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
str=[{'ANMF with KA-E'},{'ANMF with KA-PE'},{'ANMF with KA-LogE'},...
    {'ANMF with E'},{'ANMF with P'},{'ANMF with LogE'},...
    {'ANMF with CC'},{'ANMF with SFP'},{'ANMF with NMF'}];
% columnlegend(2, str);
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthWest')
% set(h_leg,'box','off')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load PdFd_4Second_s0.9_p.mat
plot(fd,circshift(Pd_ECC_mc,1),'r-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_PCC_mc,1),'g-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_LogCC_mc,1),'b-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_E_mc,1),'r-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_P_mc,1),'g-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_LogM_mc,1),'b-.','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_CC_mc,1),'k-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_SFP_mc,1),'c-*','linewidth',linewide1,'markersize',mkft)
plot(fd,circshift(Pd_R_mc,1),'m-.*','linewidth',linewide1,'markersize',mkft)

xlabel('Normalized Doppler of Target','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
str=[{'ANMF with KA-E'},{'ANMF with KA-PE'},{'ANMF with KA-LogE'},...
    {'ANMF with E'},{'ANMF with P'},{'ANMF with LogE'},...
    {'ANMF with CC'},{'ANMF with SFP'},{'ANMF with NMF'}];
% columnlegend(2, str);
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthWest')
% set(h_leg,'box','off')
