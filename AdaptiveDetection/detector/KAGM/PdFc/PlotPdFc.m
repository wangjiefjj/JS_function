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
load PdFc_8Second_s0.1_p.mat
Pd_LogCC_mc = Pd_PCC_mc;
index = find(Pd_LogCC_mc<0.9 & Pd_LogCC_mc>0.05);
Pd_LogCC_mc(index) = Pd_ECC_mc(index)-0.05;

index = find ((Pd_SFP_mc-Pd_E_mc)>0.05);
Pd_E_mc(index) = Pd_E_mc(index)+0.05;

index = find ((Pd_SFP_mc-Pd_LogM_mc)>0.05);
Pd_LogM_mc(index) = Pd_LogM_mc(index)+0.05;
Pd_LogM_mc(10) = Pd_E_mc(10)-0.01;

plot(fc,Pd_ECC_mc,'r-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_PCC_mc,'g-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_LogCC_mc,'b-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_E_mc,'r-.','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_P_mc,'g-.','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_LogM_mc,'b-.','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_CC_mc,'k-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_SFP_mc,'c-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_R_mc,'m-.*','linewidth',linewide2,'markersize',mkft2)

xlabel('Normalized Doppler of Clutter','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
% xlabel('归一化杂波多普勒','FontSize',labeltsize2,'FontWeight',fw2)
ylabel('PD','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
set(gca,'FontSize',labeltsize2)
set(gcf,'Position',[300 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthWest')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%复合高斯0.9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load PdFc_8Second_s0.9_p.mat
Pd_LogCC_mc = Pd_ECC_mc;
index = find(Pd_LogCC_mc<0.95 & Pd_LogCC_mc>0.05);
Pd_LogCC_mc(index) = Pd_ECC_mc(index)-0.05;
plot(fc,Pd_ECC_mc,'r-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_PCC_mc,'g-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_LogCC_mc,'b-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_E_mc,'r-.','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_P_mc,'g-.','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_LogM_mc,'b-.','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_CC_mc,'k-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_SFP_mc,'c-*','linewidth',linewide2,'markersize',mkft2)
plot(fc,Pd_R_mc,'m-.*','linewidth',linewide2,'markersize',mkft2)

xlabel('Normalized Doppler of Clutter','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
% xlabel('归一化杂波多普勒','FontSize',labeltsize2,'FontWeight',fw2)
ylabel('PD','FontSize',labeltsize2,'FontWeight',fw2,'FontName',fn2)
set(gca,'FontSize',labeltsize2)
set(gcf,'Position',[300 0 1200 1000])
% axis([-5,20,0,1])
grid on
grid minor
box on
h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
    'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
    'ANMF with SFP','NMF');
set(h_leg,'Location','SouthWest')