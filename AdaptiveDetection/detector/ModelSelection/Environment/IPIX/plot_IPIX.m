clc
clear
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load IPIX_recognition_Rang_19_0.mat
%%
figure()
bar(y)
set(gca, 'XTickLabel', x);
h_leg = legend('AIC','BIC','GIC (\eta=2)');
xlabel('Recognition Results')
ylabel('Percentage')
axis([0.5,3.5,0,100])
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[100 0 1200 1000])
set(h_leg,'Location','NorthWest')
grid on
box on
str=['IPIX_1','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load IPIX_recognition_Rang_19_50000.mat
%%
figure()
bar(y)
set(gca, 'XTickLabel', x);
h_leg = legend('AIC','BIC','GIC (\eta=2)');
xlabel('Recognition Results')
ylabel('Percentage')
axis([0.5,3.5,0,100])
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[100 0 1200 1000])
set(h_leg,'Location','NorthWest')
grid on
box on
str=['IPIX_2','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
