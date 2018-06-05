clc
clear
close all
FontSize = 30;
markersize =10;
linewidth = 2;
% % % % % % % % % % % % % % % ºÏ≤‚–‘ƒ‹
load ('.\Data\p_new_CCIter_ANMF_1N_s1.3.mat')
figure(1);
hold on
plot(SNRout,Pd_SCM_mc,'r','linewidth',2,'markersize',10)
plot(SNRout,Pd_CC_mc,'b','linewidth',2,'markersize',10)
plot(SNRout,Pd_ML_mc,'g','linewidth',2,'markersize',10)
plot(SNRout,Pd_CCIter_mc,'k','linewidth',2,'markersize',10)
load ('.\Data\p_new_CCIter_ANMF_2N_s1.3.mat')
plot(SNRout,Pd_SCM_mc,'r-*','linewidth',2,'markersize',10)
plot(SNRout,Pd_CC_mc,'b-*','linewidth',2,'markersize',10)
plot(SNRout,Pd_ML_mc,'g-*','linewidth',2,'markersize',10)
plot(SNRout,Pd_CCIter_mc,'k-*','linewidth',2,'markersize',10)
h_leg = legend('ANMF with SCM, K=N',...
'ANMF with CC, K=N','ANMF with ML, K=N','ANMF with KA-ICE, K=N',...
'ANMF with SCM, K=2N',...
'ANMF with CC, K=2N','ANMF with ML, K=2N','ANMF with KA-ICE, K=2N');
% xlabel('SNR/dB','FontSize',20)
xlabel({'\fontsize{30}SNR';'\fontsize{40}f'})
ylabel('PD','FontSize',FontSize)
set(gca,'FontSize',FontSize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthEast')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IPIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % ºÏ≤‚–‘ƒ‹
% load ('.\Data\R_8_new_CCIter_IPIX_19980223_170435_1N.mat')
% figure(1);
% hold on
% plot(SNRout,Pd_SCM_mc,'r','linewidth',2,'markersize',10)
% plot(SNRout,Pd_CC_mc,'b','linewidth',2,'markersize',10)
% plot(SNRout,Pd_ML_mc,'g','linewidth',2,'markersize',10)
% plot(SNRout,Pd_CCIter_mc,'k','linewidth',2,'markersize',10)
% load ('.\Data\R_8_new_CCIter_IPIX_19980223_170435_2N.mat')
% plot(SNRout,Pd_SCM_mc,'r-*','linewidth',2,'markersize',10)
% plot(SNRout,Pd_CC_mc,'b-*','linewidth',2,'markersize',10)
% plot(SNRout,Pd_ML_mc,'g-*','linewidth',2,'markersize',10)
% plot(SNRout,Pd_CCIter_mc,'k-*','linewidth',2,'markersize',10)
% h_leg = legend('ANMF with SCM, K=N',...
% 'ANMF with CC, K=N','ANMF with ML, K=N','ANMF with KA-ICE, K=N',...
% 'ANMF with SCM, K=2N',...
% 'ANMF with CC, K=2N','ANMF with ML, K=2N','ANMF with KA-ICE, K=2N');
% 
% xlabel('SNR/dB','FontSize',20)
% ylabel('Pd','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
% %%%%%%%%%%%%%%%%%%%%%%%alpha÷µ%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear
% close all
% load ('.\Data\p_CCIter_Alpha_1N.mat')
% % figure(1);
% hold on
% plot(sigma_t,m_alpha,'b','linewidth',2,'markersize',10)
% plot(sigma_t,m_alpha_ML,'g','linewidth',2,'markersize',10)
% plot(sigma_t,m_alpha_iter,'k','linewidth',2,'markersize',10)
% % clear
% load ('.\Data\p_CCIter_Alpha_2N.mat')
% plot(sigma_t,m_alpha,'b-*','linewidth',2,'markersize',10)
% plot(sigma_t,m_alpha_ML,'g-*','linewidth',2,'markersize',10)
% plot(sigma_t,m_alpha_iter,'k-*','linewidth',2,'markersize',10)
% 
% h_leg = legend('\alpha_{CC}, K=N','\alpha_{ML}, K=N','\alpha_{KA-ICE}, K=N',...
%     '\alpha_{CC}, K=2N','\alpha_{ML}, K=2N','\alpha_{KA-ICE}, K=2N');
% 
% % xlabel('\sigma^2','FontSize',20)
% % xlabel({'\fontsize{35}\sigma^2';'\fontsize{30}a'})
% xlabel({'\fontsize{35}\sigma^2';'\fontsize{30}b'})
% ylabel('\alpha','FontSize',FontSize)
% set(gca,'FontSize',FontSize)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% axis([0,10,0,1])
% grid on
% box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ŒÛ≤Ó%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear
% close all
% load ('.\Data\error_p1N.mat')
% figure(1);
% hold on
% plot(sigma_t,m_errorRSCM/100,'r','linewidth',2,'markersize',10)
% plot(sigma_t,m_errorRCC/100,'b','linewidth',2,'markersize',10)
% plot(sigma_t,m_errorRCCML/100,'g','linewidth',2,'markersize',10)
% plot(sigma_t,m_errorRCCIter/100,'k','linewidth',2,'markersize',10)
% % plot(sigma_t(1:11),m_errorRKA(1:11)/100,'y','linewidth',2,'markersize',10)
% % clear
% load ('.\Data\error_p2N.mat')
% plot(sigma_t,m_errorRSCM/100,'r-*','linewidth',2,'markersize',10)
% plot(sigma_t,m_errorRCC/100,'b-*','linewidth',2,'markersize',10)
% plot(sigma_t,m_errorRCCML/100,'g-*','linewidth',2,'markersize',10)
% plot(sigma_t,m_errorRCCIter/100,'k-*','linewidth',2,'markersize',10)
% h_leg = legend('R_{SCM}, K=N','R_{CC}, K=N','R_{ML}, K=N','R_{KA-ICE}, K=N',...
%     'R_{SCM}, K=2N','R_{CC}, K=2N','R_{ML}, K=2N','R_{KA-ICE}, K=2N');
% % xlabel({'\fontsize{35}\sigma^2';'\fontsize{30}a'})
% xlabel({'\fontsize{35}\sigma^2';'\fontsize{30}b'})
% ylabel('Error_{NFN}','FontSize',FontSize)
% set(gca,'FontSize',FontSize)
% set(gcf,'Position',[700 0 1200 1000])
% set(h_leg,'Location','SouthEast')
% grid on
% box on
%%%%%%%%%%%%IPIXæ‡¿Î∂‡∆’¿’Õº%%%%%%%%%%%%
% clc
% clear
% close all
% load ('.\Data\19980223_170435_IPIX.mat')
% [M,L] = size(sig);
% LL = 1:L;
% MM = linspace(-0.5,0.5,M);
% [X,Y]=meshgrid(LL,MM);
% MTD = abs(fftshift(fft(sig,[],1)));
% [x,y] = max(MTD);
% 
% figure()
% plot3(8,MM(y(8)),1,'ro','markersize',10,'MarkerFaceColor','r')
% str = [' Range cell: %.f \n Normalized Doppler: %.6f \n Normalized amplitude: %.3f'];
% text(8,MM(y(8)),1,sprintf(str,8,MM(y(8)),1),...
%     'VerticalAlignment','bottom','FontSize',16)
% hold on
% mesh(X,Y,MTD/max(max(MTD)));
% ylabel('Normalized Doppler','FontSize',20)
% xlabel('Range cell','FontSize',20)
% zlabel('Normalized amplitude','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])
% grid on
% box on
% axis([1,34,-0.5,0.5,0,1])


% figure()
% imagesc(1:L,linspace(-0.5,0.5,M),MTD/max(max(MTD)));
% ylabel('Normalized Doppler','FontSize',20)
% xlabel('Range cell','FontSize',20)
% zlabel('Normalized amplitude','FontSize',20)
% set(gca,'FontSize',20)
% set(gcf,'Position',[700 0 1200 1000])