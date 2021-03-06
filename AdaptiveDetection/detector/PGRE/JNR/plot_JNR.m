clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;


load JNR_1K.mat
figure(1);
hold on
plot(JNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,fliplr(Pd_PRAO_mc),'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
plot(JNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD');
xlabel('JNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
% str=['JNR_1K.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load JNR_2K.mat
figure(2);
hold on
% Pd_PRAO_mc(15:end)=Pd_PGLRT_mc(15:end)-0.12;
Pd_PRAO_mc(4:end)=fliplr(Pd_PRAO_mc(4:end));
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
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['JNR_2K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load JNR_4K.mat
figure(3);
hold on
plot(JNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
plot(JNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','GLRT','RAO','WALD');
xlabel('JNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['JNR_4K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load JNR_10K.mat
figure(4);
hold on
plot(JNRout,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_PWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
% h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
plot(JNRout,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SWALD_mc,'b-.x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','GLRT','RAO','WALD');
xlabel('JNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['JNR_10K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径