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
load Confusion_L20.mat
%% 
%% 综合来看ABIC 和 GIC rho=2  比较好 
%%Classification results for H1 hypothesis assuming
%%真实环境为H1时判断为：H1，H2，H3
y_H1 =[Confusion_AIC(1,1),Confusion_AICc(1,1),Confusion_ABIC(1,1),...
    Confusion_GIC1(1,1),Confusion_GIC2(1,1)];
y_H2 =[Confusion_AIC(1,2),Confusion_AICc(1,2),Confusion_ABIC(1,2),...
    Confusion_GIC1(1,2),Confusion_GIC2(1,2)];
y_H3 =[Confusion_AIC(1,3),Confusion_AICc(1,3),Confusion_ABIC(1,3),...
    Confusion_GIC1(1,3),Confusion_GIC2(1,3)];
y = [y_H1;y_H2;y_H3];
% y = y';
x={'H1','H2','H3'};
%%
figure()
bar(y)
set(gca, 'XTickLabel', x);
h_leg = legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('Classification results')
ylabel('Percentage')
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[100 0 1200 1000])
grid on
box on
% str=['ConfusionL12_H1.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%Classification results for H2 hypothesis assuming
%%真实环境为H2时判断为：H1，H2，H3
y_H1 =[Confusion_AIC(2,1),Confusion_AICc(2,1),Confusion_ABIC(2,1),...
    Confusion_GIC1(2,1),Confusion_GIC2(2,1)];
y_H2 =[Confusion_AIC(2,2),Confusion_AICc(2,2),Confusion_ABIC(2,2),...
    Confusion_GIC1(2,2),Confusion_GIC2(2,2)];
y_H3 =[Confusion_AIC(2,3),Confusion_AICc(2,3),Confusion_ABIC(2,3),...
    Confusion_GIC1(2,3),Confusion_GIC2(2,3)];
y = [y_H1;y_H2;y_H3];
% y = y';
x={'H1','H2','H3'};
%%
figure()
bar(y)
set(gca, 'XTickLabel', x);
h_leg = legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('Classification results')
ylabel('Percentage')
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[300 0 1200 1000])
grid on
box on
% str=['ConfusionL12_H2.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%
%%Classification results for H3 hypothesis assuming
%%真实环境为H3时判断为：H1，H2，H3
y_H1 =[Confusion_AIC(3,1),Confusion_AICc(3,1),Confusion_ABIC(3,1),...
    Confusion_GIC1(2,1),Confusion_GIC2(2,1)];
y_H2 =[Confusion_AIC(3,2),Confusion_AICc(3,2),Confusion_ABIC(3,2),...
    Confusion_GIC1(3,2),Confusion_GIC2(3,2)];
y_H3 =[Confusion_AIC(3,3),Confusion_AICc(3,3),Confusion_ABIC(3,3),...
    Confusion_GIC1(3,3),Confusion_GIC2(3,3)];
y = [y_H1;y_H2;y_H3];
% y = y';
x={'H1','H2','H3'};
%%
figure()
bar(y)
set(gca, 'XTickLabel', x);
h_leg = legend('AIC','AICc','ABIC','GIC(\rho=2)','GIC(\rho=4)');
xlabel('Classification results')
ylabel('Percentage')
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(h_leg,'Location','NorthWest')
set(gcf,'Position',[500 0 1200 1000])
grid on
box on
% str=['ConfusionL12_H3.eps'];
% print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%bar3，似乎不好使
%AIC
% figure()
% bar3(Confusion_ABIC)

