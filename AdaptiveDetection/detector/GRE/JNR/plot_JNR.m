clc
clear 
close all
labeltsize=25;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=5;
mkft = 15;
load JNR.mat
figure(1);
hold on
plot(JNRout,Pd_SGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',5)
plot(JNRout,Pd_SRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(JNRout,Pd_SWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('GLRT','RAO','WALD');
xlabel('JNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([0,42,0,1])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthWest')
grid on
box on
str=['JNR.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径