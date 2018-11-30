clc
clear 
close all
labeltsize=25;
fw = 'normal'; %%�Ƿ�Ӵ�б��֮��
fn='Times New Roman';
linewide1=5;
mkft = 15;
load fc.mat
figure(1);
hold on
plot(-fc,Pd_SGLRT_mc,'k-o','linewidth',linewide1,'MarkerSize',5)
plot(-fc,Pd_SRAO_mc,'r-x','linewidth',linewide1,'MarkerSize',mkft)
plot(-fc,Pd_SWALD_mc,'b-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg=legend('GLRT','RAO','WALD');
xlabel('f_{c}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
axis([-0.5,0.5,0,1])
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','NorthWest')
grid on
box on
str=['fc.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��