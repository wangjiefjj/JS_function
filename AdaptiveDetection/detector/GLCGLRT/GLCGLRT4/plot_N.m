clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下数据量的关系%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load Pd_CLGLRT4_2_1Kmu1lambda3s0.1o1_p.mat
% Pd_CLGLRT_mc(1:4)=Pd_CLGLRT_mc(1:4)*0.8;
% Pd_CLGLRT_mc(5:end)=Pd_CLGLRT_mc(1:end-4);
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_1Kmu1lambda3_p.mat
plot(SNRout,Pd_SAML_mc,'r-p','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_p1.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
    
figure(2)
hold on
load Pd_CLGLRT4_2_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_2Kmu1lambda3_p.mat
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_p2.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%高斯环境下%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
load Pd_CLGLRT4_2_1Kmu1lambda3s0.1o1_g.mat
% Pd_CLGLRT_mc(1:4)=Pd_CLGLRT_mc(1:4)*0.8;
% Pd_CLGLRT_mc(5:end)=Pd_CLGLRT_mc(1:end-4);
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_1Kmu1lambda3_g.mat
plot(SNRout,Pd_SAML_mc,'r-p','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_g1.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径


figure(4)
hold on
load Pd_CLGLRT4_2_2Kmu1lambda3s0.1o1_g.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_2Kmu1lambda3_g.mat
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_g2.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
% pause(2)
% close all

figure(5)
hold on
load Pd_CLGLRT4_2_3Kmu1lambda3s0.1o1_g.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_3Kmu1lambda3_g.mat
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_g3.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径

figure(5)
hold on
load Pd_CLGLRT4_2_10Kmu1lambda3s0.1o1_g.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);%Pd_CLGLRT_mc/2+Pd_KGLRTCC_mc/2
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_g4.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径