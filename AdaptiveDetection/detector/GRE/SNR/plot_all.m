clc
clear 
close all
labeltsize=25;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=5;
mkft = 15;
load 2N.mat
figure(1);
hold on
plot(SNRout,Pd_SGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',5)
plot(SNRout,Pd_SRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('GLRT','RAO','WALD');
xlabel('SNR(dB)','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1300 1000])
set(h_leg,'Location','NorthWest')
grid on
box on
str=['SNR2N.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 10N.mat
figure(2);
hold on
plot(SNRout,Pd_SGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',5)
plot(SNRout,Pd_SRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(SNRout,Pd_SWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('GLRT','RAO','WALD');
xlabel('SNR(dB)','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1300 1000])
set(h_leg,'Location','NorthWest')
grid on
box on
str=['SNR10N.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
pause(1)
close all