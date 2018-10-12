clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mis_1K35dB.mat
figure(1);
hold on
plot(cos2,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_PWALD_mc,'b-*','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_SWALD_mc,'b-.*','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD');
xlabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 500])
set(h_leg,'Location','NorthWest')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mis_2K35dB.mat
figure(2);
hold on
plot(cos2,Pd_PGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_PRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_PWALD_mc,'b-*','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_SGLRT_mc,'k-.o','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_SRAO_mc,'r-.x','linewidth',linewide1,'MarkerSize',mkft)
plot(cos2,Pd_SWALD_mc,'b-.*','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD','SGLRT','SRAO','SWALD');
xlabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 500])
set(h_leg,'Location','NorthWest')
grid on
box on