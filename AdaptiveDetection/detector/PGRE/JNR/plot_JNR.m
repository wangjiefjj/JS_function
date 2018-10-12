clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
load JNR_1K.mat
figure(1);
hold on
plot(JNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
plot(JNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD');
xlabel('JNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 500])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load JNR_2K.mat
figure(2);
hold on
plot(JNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
plot(JNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD');
xlabel('JNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 500])
set(h_leg,'Location','SouthEast')
grid on
box on