clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%非高斯环境下，RKA=单位阵%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
hold on
SNRout_real=0:1:25; % 输出SNR
L = length(SNRout_real);
%%1N 
load Pd_CLGLRT_revision_one_1Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_CLGLRT_mc(5:L+4),'k-p','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-p','linewidth',2,'markersize',10);
load Pd_CLGLRT4_2_1Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b-p','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c-p','linewidth',2,'markersize',10);
load Pd_CLGLRT_revision_saml_1Kmu1lambda3_p.mat
plot(SNRout_real,Pd_SAML_mc(1:L),'r-p','linewidth',2,'markersize',10);
%%2N
load Pd_CLGLRT_revision_one_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_CLGLRT_mc(2:L+1),'k->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g->','linewidth',2,'markersize',10);
load Pd_CLGLRT4_2_2Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b->','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c->','linewidth',2,'markersize',10);
load Pd_CLGLRT_revision_saml_2Kmu1lambda3_p.mat
plot(SNRout_real,Pd_SAML_mc(1:L),'r->','linewidth',2,'markersize',10);
%%3N
load Pd_CLGLRT_revision_one_3Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_CLGLRT_mc(2:L+1),'k-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTCC_mc(1:L),'g-o','linewidth',2,'markersize',10);
load Pd_CLGLRT4_2_3Kmu1lambda3s0.1o1_p.mat
plot(SNRout_real,Pd_KGLRT_mc(1:L),'b-o','linewidth',2,'markersize',10);
plot(SNRout_real,Pd_KGLRTNSCM_mc(1:L),'c-o','linewidth',2,'markersize',10);
load Pd_CLGLRT_revision_saml_3Kmu1lambda3_p.mat
plot(SNRout_real,Pd_SAML_mc(1:L),'r-o','linewidth',2,'markersize',10);
h_leg = legend('GLC-GLRT,K=N','1S-GLRT with CC,K=N',...
               '1S-GLRT with SCM,K=N','1S-GLRT with NSCM,K=N','1S-GLRT with SFPE,K=N',...
               'GLC-GLRT,K=2N','1S-GLRT with CC,K=2N',...
               '1S-GLRT with SCM,K=2N','1S-GLRT with NSCM,K=2N','1S-GLRT with SFPE,K=2N',...
                'GLC-GLRT,K=3N','1S-GLRT with CC,K=3N',...
               '1S-GLRT with SCM,K=3N','1S-GLRT with NSCM,K=3N','1S-GLRT with SFPE,K=3N');
xlabel('SNR/dB','FontSize',20)%,'FontWeight',fw,'FontName',fn
ylabel('Pd','FontSize',20)%,'FontWeight',fw,'FontName',fn
set(gca,'FontSize',20) %,'FontWeight',fw,'FontName',fn
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on