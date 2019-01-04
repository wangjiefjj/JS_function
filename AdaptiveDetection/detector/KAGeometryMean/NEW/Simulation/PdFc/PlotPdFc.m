clc
clear
close all
load PdFc2_4Second_s0.1_p.mat
figure();
hold on
plot(fc,Pd_R_mc,'c.-','linewidth',1)
plot(fc,Pd_CC_mc,'b-*','linewidth',1)
plot(fc,Pd_E_mc,'k-.*','linewidth',1)
plot(fc,Pd_ECC_mc,'k-*','linewidth',1)
plot(fc,Pd_LogM_mc,'r.-','linewidth',1)
plot(fc,Pd_LogCC_mc,'r-*','linewidth',1)
plot(fc,Pd_P_mc,'g.-','linewidth',1)
load PdFcPCC_4Second_s0.1_p.mat
plot(fc,Pd_PCC_mc,'g-*','linewidth',1)
load PdFc2_4Second_s0.1_p.mat
plot(fc,Pd_SFP_mc,'c-*','linewidth',1)
h_leg = legend('NMF','ANMF with CC',...
    'ANMF with E','ANMF with ECC','ANMF with LogM','ANMF with LogCC',...
    'ANMF with P','ANMF with PCC','SFP');