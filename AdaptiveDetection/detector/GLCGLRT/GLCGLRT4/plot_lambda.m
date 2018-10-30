clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%�Ƿ�Ӵ�б��֮��
fn='Times New Roman';
linewide1=3;
mkft = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%��ͬlambda�µļ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
load Pd_CLGLRT4_2_2Kmu1lambda1s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-p','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-p','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_2Kmu1lambda1_p.mat
plot(SNRout,Pd_SAML_mc,'r-p','linewidth',linewide1,'markersize',mkft);
h_leg = legend...
       ('GLC-GLRT',        '1S-GLRT with CC ',...
       '1S-GLRT with SCM ','1S-GLRT with NSCM ',...
       '1S-GLRT with SFPE ');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on 
box on
str=['lambda1.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��

figure(2)
hold on
load Pd_CLGLRT4_2_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_2Kmu1lambda3_p.mat
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);
h_leg = legend...
       ('GLC-GLRT',        '1S-GLRT with CC ',...
       '1S-GLRT with SCM ','1S-GLRT with NSCM ',...
       '1S-GLRT with SFPE ');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on 
box on
str=['lambda3.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��

figure(3)
hold on
load Pd_CLGLRT4_2_2Kmu1lambda4s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b->','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c->','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_2Kmu1lambda3_p.mat
plot(SNRout,Pd_SAML_mc,'r->','linewidth',linewide1,'markersize',mkft);
h_leg = legend...
       ('GLC-GLRT',        '1S-GLRT with CC ',...
       '1S-GLRT with SCM ','1S-GLRT with NSCM ',...
       '1S-GLRT with SFPE ');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on 
box on
str=['lambda4.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��

figure(4)
hold on
load Pd_CLGLRT4_2_2Kmu1lambda5s0.1o1_p.mat
plot(SNRout,Pd_CLGLRT_mc,'k-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTCC_mc,'g-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRT_mc,'b-o','linewidth',linewide1,'markersize',mkft);
plot(SNRout,Pd_KGLRTNSCM_mc,'c-o','linewidth',linewide1,'markersize',mkft);
load Pd_CLGLRT_revision_saml_2Kmu1lambda5_p.mat
plot(SNRout,Pd_SAML_mc,'r-o','linewidth',linewide1,'markersize',mkft);
h_leg = legend...
       ('GLC-GLRT',        '1S-GLRT with CC ',...
       '1S-GLRT with SCM ','1S-GLRT with NSCM ',...
       '1S-GLRT with SFPE ');
xlabel('SCR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on 
box on
str=['lambda5.eps'];
print(gcf,'-depsc',str)   %����Ϊpng��ʽ��ͼƬ����ǰ·��