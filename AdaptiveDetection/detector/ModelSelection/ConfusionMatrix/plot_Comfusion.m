clc
clear
close all
labeltsize=35;
fw = 'normal'; %%�Ƿ�Ӵ�б��֮��
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load Confusion_L20.mat
%% 
%% �ۺ�����ABIC �� GIC rho=2  �ȽϺ� 
%%Classification results for H1 hypothesis assuming
%%��ʵ����ΪH1ʱ�ж�Ϊ��H1��H2��H3
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
% print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%Classification results for H2 hypothesis assuming
%%��ʵ����ΪH2ʱ�ж�Ϊ��H1��H2��H3
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
% print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%
%%Classification results for H3 hypothesis assuming
%%��ʵ����ΪH3ʱ�ж�Ϊ��H1��H2��H3
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
% print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%bar3���ƺ�����ʹ
%AIC
% figure()
% bar3(Confusion_ABIC)

