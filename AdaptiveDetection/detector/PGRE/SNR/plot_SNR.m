clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
load 1K.mat
figure(1);
hold on
plot(SNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['Pd_1K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 2K.mat
figure(2);
hold on
plot(SNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('PGLRT','PRAO','PWALD','SGLRT','SRAO','SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['Pd_2K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load 10K.mat
% figure(3);
% hold on
% plot(SNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
% plot(SNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('PGLRT','PRAO','PWALD','SGLRT','SRAO','SWALD');
% xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% set(gca,'FontSize',labeltsize)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on