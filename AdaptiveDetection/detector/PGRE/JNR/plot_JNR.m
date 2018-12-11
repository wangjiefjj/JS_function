clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%�Ƿ�Ӵ�б��֮��
fn='Times New Roman';
linewide1=3;
mkft = 15;


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
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['JNR_1K.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
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
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['JNR_2K.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��