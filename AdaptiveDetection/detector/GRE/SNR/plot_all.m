clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
load 2K.mat
figure(1);
hold on
plot(SNRout,Pd_SGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('SGLRT','SRAO','SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 10K.mat
figure(2);
hold on
plot(SNRout,Pd_SGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('SGLRT','SRAO','SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('Pd','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on