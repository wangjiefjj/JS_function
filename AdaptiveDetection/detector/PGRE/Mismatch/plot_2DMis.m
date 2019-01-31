clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mis_1K35dB.mat
figure(1);
hold on
plot(1-cos2,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PWALD_mc,'b-*','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SWALD_mc,'b-.*','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','GLRT','RAO','WALD');
xlabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','NorthEast')
grid on
box on
str=['Mis2D_1K35dB.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mis_2K35dB.mat
figure(2);
hold on
plot(1-cos2,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PWALD_mc,'b-*','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SWALD_mc,'b-.*','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','GLRT','RAO','WALD');
xlabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','NorthEast')
grid on
box on
str=['Mis2D_2K35dB.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mis_4K35dB.mat
figure(3);
hold on
plot(1-cos2,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PWALD_mc,'b-*','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SWALD_mc,'b-.*','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','GLRT','RAO','WALD');
xlabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthWest')
grid on
box on
str=['Mis2D_4K35dB.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mis_10K35dB.mat
figure(4);
hold on
plot(1-cos2,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_PWALD_mc,'b-*','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(1-cos2,Pd_SWALD_mc,'b-.*','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','GLRT','RAO','WALD');
xlabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthWest')
grid on
box on
str=['Mis2D_10K35dB.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径