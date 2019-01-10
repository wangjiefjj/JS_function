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

figure(1)
hold on
load IPIX_19980204_224024_1N.mat
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(SNRout),max(SNRout),0,1])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['IPIX_PD1N.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
load IPIX_19980204_224024_2N.mat
plot(SNRout,Pd_SCM_mc,'k-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_CC_mc,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ML_mc,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCP_mc,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCS_mc,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_ECCT_mc,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('ANMF with SCM','ANMF with CC','ANMF with ML',...
    'ANMF with KA-P','ANMF with KA-S','ANMF with KA-T');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([min(SNRout),max(SNRout),0,1])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['IPIX_PD2N.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径