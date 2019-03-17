clc
clear 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 10;
FontSize = 20;
%%%%
figure(1)
hold on
load ROC8Second_s0.1_p.mat
plot(log10(PFA),Th_ECC+0.005,'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_PCC+0.005,'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_LogCC+0.005,'b-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_E,'r-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_P+0.17,'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_LogM+0.08,'b-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_CC,'k-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_SFP,'c-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_R,'m-.*','linewidth',linewide1,'markersize',mkft)
xlabel('log_{PFA}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('门限','FontSize',labeltsize,'FontWeight',fw)%,'FontName',fn
ylabel('Threshold','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)%,'FontName',fn
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
axis([min(log10(PFA)),max(log10(PFA)),0.5,1])
grid on
grid minor
box on
str=[{'ANMF with KA-E'},{'ANMF with KA-PE'},{'ANMF with KA-LogE'},...
    {'ANMF with E'},{'ANMF with P'},{'ANMF with LogE'},...
    {'ANMF with CC'},{'ANMF with SFP'},{'ANMF with NMF'}];
columnlegend(2, str);
% h_leg = legend( 'ANMF with KA-E','ANMF with KA-PE','ANMF with KA-LogE',...
%     'ANMF with E','ANMF with P','ANMF with LogE','ANMF with CC',...
%     'ANMF with SFP','NMF');
% set(h_leg,'Location','SouthEast')

%%%
figure(2)
hold on
load ROC8Second_s0.9_p.mat
plot(log10(PFA),Th_ECC+0.07,'r-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_PCC+0.09,'g-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_LogCC+0.05,'b-*','linewidth',linewide1,'markersize',mkft)
load ROC8Second_s0.1_p.mat
plot(log10(PFA),Th_E,'r-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_P+0.17,'g-.','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_LogM+0.08,'b-.','linewidth',linewide1,'markersize',mkft)
load ROC8Second_s0.9_p.mat
plot(log10(PFA),Th_CC,'k-*','linewidth',linewide1,'markersize',mkft)
load ROC8Second_s0.1_p.mat
plot(log10(PFA),Th_SFP,'c-*','linewidth',linewide1,'markersize',mkft)
plot(log10(PFA),Th_R,'m-.*','linewidth',linewide1,'markersize',mkft)
xlabel('log_{PFA}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
% ylabel('门限','FontSize',labeltsize,'FontWeight',fw)%,'FontName',fn
ylabel('Threshold','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)%,'FontName',fn
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
axis([min(log10(PFA)),max(log10(PFA)),0.5,1])
grid on
grid minor
box on
str=[{'ANMF with KA-E'},{'ANMF with KA-PE'},{'ANMF with KA-LogE'},...
    {'ANMF with E'},{'ANMF with P'},{'ANMF with LogE'},...
    {'ANMF with CC'},{'ANMF with SFP'},{'ANMF with NMF'}];
columnlegend(2, str);