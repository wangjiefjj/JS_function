clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%实测数据IPIX%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
FontSize = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% hold on
% load Pd_10Second_19980223_170435_IPIX.mat
% plot(SNRout,Pd_CC_mc,'r-','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_ML_mc,'r-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_NSCM_mc,'g-','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_ECC_mc,'g-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_LogM_mc,'b-','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_P_mc,'k-','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_PCC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
% load Pd_16Second_19980223_170435_IPIX.mat
% plot(SNRout,Pd_CC_mc,'r-.','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_ML_mc,'r-.*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_NSCM_mc,'g-.','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_ECC_mc,'g-.*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_LogCC_mc,'b-.*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_P_mc,'k-.','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_PCC_mc,'k-.*','linewidth',linewide1,'markersize',mkft)
% 
% 
% 
% xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[700 0 1200 1000])
% % axis([-5,15,0,1])
% grid on
% box on
% h_leg = legend('ANMF with CC, K=10','ANMF with ML, K=10','ANMF with NSCM, K=10',...
%     'ANMF with ECC, K=10','ANMF with LogM, K=10','ANMF with LogCC, K=10',...
%     'ANMF with P, K=10','ANMF with PCC, K=10','ANMF with CC, K=16',...
%     'ANMF with ML, K=16','ANMF with NSCM, K=16',...
%     'ANMF with ECC, K=16','ANMF with LogM, K=16','ANMF with LogCC, K=16',...
%     'ANMF with P, K=16','ANMF with PCC, K=16');
% set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
load Pd_10Second19980204_224024_IPIX.mat
plot(SNRout,Pd_CC_mc,'r-','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_ML_mc,'r-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_NSCM_mc,'g-','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_ECC_mc,'g-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogM_mc,'b-','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_P_mc,'k-','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_PCC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
load Pd_16Second19980204_224024_IPIX.mat
plot(SNRout,Pd_CC_mc,'r-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_ML_mc,'r-.*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_NSCM_mc,'g-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_ECC_mc,'g-.*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_LogCC_mc,'b-.*','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_P_mc,'k-.','linewidth',linewide1,'markersize',mkft)
plot(SNRout,Pd_PCC_mc,'k-.*','linewidth',linewide1,'markersize',mkft)



xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,15,0,1])
grid on
box on
h_leg = legend('ANMF with CC, K=10','ANMF with ML, K=10','ANMF with NSCM, K=10',...
    'ANMF with ECC, K=10','ANMF with LogM, K=10','ANMF with LogCC, K=10',...
    'ANMF with P, K=10','ANMF with PCC, K=10','ANMF with CC, K=16',...
    'ANMF with ML, K=16','ANMF with NSCM, K=16',...
    'ANMF with ECC, K=16','ANMF with LogM, K=16','ANMF with LogCC, K=16',...
    'ANMF with P, K=16','ANMF with PCC, K=16');
set(h_leg,'Location','SouthEast')
%%%%%%%%%%%%IPIX距离多普勒图%%%%%%%%%%%%
clc
clear
close all
load ('19980204_224024_IPIX.mat')%19980204_224024_IPIX%19980223_170435_IPIX
[M,L] = size(sig);
LL = 1:L;
MM = linspace(-0.5,0.5,M);
[X,Y]=meshgrid(LL,MM);
MTD = abs(fftshift(fft(sig,[],1)));
[x,y] = max(MTD);
[~,max_Range] = max(x);
figure()
plot3(max_Range,MM(y(max_Range)),1,'ro','markersize',10,'MarkerFaceColor','r')
str = [' Range cell: %.f \n Normalized Doppler: %.6f \n Normalized amplitude: %.3f'];
text(max_Range,MM(y(max_Range)),1,sprintf(str,max_Range,MM(y(max_Range)),1),...
    'VerticalAlignment','bottom','FontSize',25)
hold on
mesh(X,Y,MTD/max(max(MTD)));
ylabel('Normalized  Doppler','FontSize',FontSize)
xlabel('Range cell','FontSize',FontSize)
zlabel('Normalized amplitude','FontSize',FontSize)
set(gca,'FontSize',FontSize)
set(gcf,'Position',[700 0 1200 1000])
grid on
box on
axis([1,34,-0.5,0.5,0,1])