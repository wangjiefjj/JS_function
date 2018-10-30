clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下，RKA=单位阵%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
SNRout_real=0:1:25; % 输出SNR
L = length(SNRout_real);
figure(1)
hold on
%%1N 
load Pd_CLGLRT_revision_one_1Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_CLGLRT_mc(5:L+4),'k-p','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-p','linewidth',2,'markersize',10);
load Pd_CLGLRT4_2_1Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b-p','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c-p','linewidth',2,'markersize',10);
load Pd_CLGLRT_revision_saml_1Kmu1lambda3_p.mat
plot(SNRout_real,Pd_SAML_mc(1:L),'r-p','linewidth',2,'markersize',10);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)%,'FontWeight',fw,'FontName',fn
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)%,'FontWeight',fw,'FontName',fn
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn) %,'FontWeight',fw,'FontName',fn
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_one1.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径

figure(2)
hold on
%%2N
load Pd_CLGLRT_revision_one_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_CLGLRT_mc(2:L+1),'k->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g->','linewidth',2,'markersize',10);
load Pd_CLGLRT4_2_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c->','linewidth',2,'markersize',10);
load Pd_CLGLRT_revision_saml_2Kmu1lambda3_p.mat
plot(SNRout_real,Pd_SAML_mc(1:L),'r->','linewidth',2,'markersize',10);
h_leg = legend('GLC-GLRT','1S-GLRT with CC',...
               '1S-GLRT with SCM','1S-GLRT with NSCM','1S-GLRT with SFPE');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)%,'FontWeight',fw,'FontName',fn
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)%,'FontWeight',fw,'FontName',fn
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn) %,'FontWeight',fw,'FontName',fn
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
str=['PD_one2.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径