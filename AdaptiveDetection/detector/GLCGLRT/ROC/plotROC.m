clc
clear 
close all
load Pd_CLGLRT2_ROC22Kmu1lambda3s0.1o1_p.mat
load Pd_ROC_aml2Kmu1lambda3s0.1o1_p.mat
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;

figure
hold on
%%5dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,1),'k-s','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,1),'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRT_Mlti_mc(:,1),'b->','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,1),'c-*','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_SAML_Mlti_mc(:,1),'r-*','linewidth',linewide1,'MarkerSize',mkft)

%%10dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,2),'k-s','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,2),'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRT_Mlti_mc(:,2),'b->','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,2),'c-*','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_SAML_Mlti_mc(:,2),'r-*','linewidth',linewide1,'MarkerSize',mkft)
%%15dB
plot(PFA,Pd_CLGLRT_Mlti_mc(:,3),'k-s','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTCC_Mlti_mc(:,3),'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRT_Mlti_mc(:,3),'b->','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_KGLRTNSCM_Mlti_mc(:,3),'c-*','linewidth',linewide1,'MarkerSize',mkft)
plot(PFA,Pd_SAML_Mlti_mc(:,3),'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('GLC-GLRT, K=2N','1S-GLRT with CC, K=2N','1S-GLRT with SCM, K=2N',...
               '1S-GLRT with NSCM, K=2N','1S-GLRT with SFPE, K=2N');
xlabel('PFA','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('PD','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','SouthEast','FontWeight',fw,'FontName',fn)
set(gca,'Xscale','log')
axis([1e-4,1e-1,0,1])
grid on
box on