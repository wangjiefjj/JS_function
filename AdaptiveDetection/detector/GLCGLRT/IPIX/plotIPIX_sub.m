clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load Pd_CLGLRT_2_19980223_170435_IPI_1Ks0.mat
% Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_CLGLRT_mc,'k-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-*','linewidth',linewide1,'markersize',mkft);
% load Pd_CLGLRT_2_19980223_170435_IPI_1Ks0.mat
% Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
load Pd_CLGLRT_one_19980223_170435_IPI_1Ks0.mat
plot(SNRout,Pd_KGLRT_mc,'b-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_SAML_mc,'r-*','linewidth',linewide1,'markersize',mkft);
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,15,0,1])
grid on
box on
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
set(h_leg,'Location','SouthEast')
str=['IPIX_PD1.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径

figure(2)
hold on
load Pd_CLGLRT_2_19980223_170435_IPI_2Ks0.mat
% Pd_CLGLRT_mc(1:end-1)=Pd_CLGLRT_mc(2:end);
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_one_19980223_170435_IPI_2Ks0.mat
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);

xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,15,0,1])
grid on
box on
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
set(h_leg,'Location','SouthEast')
str=['IPIX_PD2.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%IPIX_one%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
hold on
load Pd_CLGLRT_one_19980223_170435_IPI_1Ks0.mat
plot(SNRout,Pd_CLGLRT_mc,'k-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-*','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_SAML_mc,'r-*','linewidth',linewide1,'markersize',mkft);
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,15,0,1])
grid on
box on
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
set(h_leg,'Location','SouthEast')
str=['IPIX_one1.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径

figure(4)
hold on
load Pd_CLGLRT_one_19980223_170435_IPI_2Ks0.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
% axis([-5,15,0,1])
grid on
box on
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
set(h_leg,'Location','SouthEast')
str=['IPIX_one2.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径