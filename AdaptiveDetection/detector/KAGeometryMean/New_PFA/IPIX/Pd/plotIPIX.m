clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%实测数据IPIX%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 10;
FontSize = 35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%19980223_170435%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% hold on
% load PD2_8Second19980223_170435_IPIX.mat
% plot(SNRout,Pd_ECC_mc,'r-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_PCC_mc,'g-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_LogCC_mc,'b-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_E_mc,'r-.','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_P_mc,'g-.','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_LogM_mc,'b-.','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_CC_mc,'k-*','linewidth',linewide1,'markersize',mkft)
% plot(SNRout,Pd_SFP_mc,'c-*','linewidth',linewide1,'markersize',mkft)
% % plot(SNRout,Pd_NSCM_mc,'m-.','linewidth',linewide1,'markersize',mkft)
% 
% xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[700 0 1200 1000])
% axis([-5,25,0,1])
% grid on
% grid minor
% box on
% h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
%     'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
%     'ANMF with SFP');
% set(h_leg,'Location','SouthEast')


%%%%%%%%%%%%IPIX距离多普勒图%%%%%%%%%%%%
load ('19980223_170435_IPIX.mat')%19980204_224024_IPIX%19980223_170435_IPIX
[M,L] = size(sig);
LL = 1:L;
MM = linspace(-0.5,0.5,M);
[X,Y]=meshgrid(LL,MM);
MTD = abs(fftshift(fft(sig,[],1),1));
MTD=MTD/max(max(MTD));
[x,y] = max(MTD);
[~,max_Range] = max(x);
figure(3)
plot3(max_Range,MM(y(max_Range)),max(max(MTD)),'ro','markersize',10,'MarkerFaceColor','r')
str = [' Range cell: %.f \n Normalized Doppler: %.6f \n Normalized amplitude: %.3f'];
text(max_Range,MM(y(max_Range)),max(max(MTD)),sprintf(str,max_Range,MM(y(max_Range)),max(max(MTD))),...
    'VerticalAlignment','bottom','FontSize',35)
hold on

mesh(X,Y,MTD);
ylabel('Normalized  Doppler','FontSize',FontSize)
xlabel('Range cell','FontSize',FontSize)
zlabel('Normalized amplitude','FontSize',FontSize)
set(gca,'FontSize',FontSize)
set(gcf,'Position',[700 0 1200 1000])
grid on
box on
axis([1,34,-0.5,0.5,0,1])