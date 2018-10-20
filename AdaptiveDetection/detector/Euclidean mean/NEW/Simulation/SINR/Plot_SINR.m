clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
load p_SINR_1N_s0.1.mat
figure(1);
hold on
plot(ft,R_PSD,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_SCM,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_CC,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCT,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCS,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCP,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('R','SCM','CC','ML','ECCT','ECCS','ECCP');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([-0.5,0.5,-20,1]);
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_SINR_1N_s0.9.mat
figure(2);
hold on
plot(ft,R_PSD,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_SCM,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_CC,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCT,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCS,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCP,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('R','SCM','CC','ML','ECCT','ECCS','ECCP');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([-0.5,0.5,-20,1]);
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_SINR_2N_s0.1.mat
figure(3);
hold on
plot(ft,R_PSD,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_SCM,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_CC,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCT,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCS,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCP,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('R','SCM','CC','ML','ECCT','ECCS','ECCP');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([-0.5,0.5,-20,1]);
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_SINR_2N_s0.9.mat
figure(4);
hold on
plot(ft,R_PSD,'r','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_SCM,'r-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_CC,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCT,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCS,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(ft,mean_SINR_ECCP,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('R','SCM','CC','ML','ECCT','ECCS','ECCP');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
axis([-0.5,0.5,-20,1]);
grid on
box on