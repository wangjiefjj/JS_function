clc
clear 
close all
labeltsize=20;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
load p_Alpha_1N.mat
figure(1);
hold on
plot(sigma_t,m_alpha,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_P,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_S,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_T,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','ECCP','ECCS','ECCT');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthWest')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_Alpha_2N.mat
figure(2);
hold on
plot(sigma_t,m_alpha,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_P,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_S,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_T,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','ECCP','ECCS','ECCT');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthWest')
grid on
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load g_Alpha_1N.mat
figure(3);
hold on
plot(sigma_t,m_alpha,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_P,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_S,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_T,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','ECCP','ECCS','ECCT');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthWest')
grid on
box on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load g_Alpha_2N.mat
figure(4);
hold on
plot(sigma_t,m_alpha,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_ML,'g-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_P,'g-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_S,'b-x','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t,m_alpha_T,'r-x','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','ECCP','ECCS','ECCT');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[700 0 1200 1000])
set(h_leg,'Location','SouthWest')
grid on
box on