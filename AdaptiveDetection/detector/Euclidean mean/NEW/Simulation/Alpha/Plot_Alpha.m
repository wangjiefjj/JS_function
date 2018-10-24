clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=2;
mkft = 10;
load p_Alpha_1N.mat
figure(1);
hold on
plot(sigma_t.^2,m_alpha,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_ML,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_P,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_S,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_T,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','KA-P','KA-S','KA-T');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','NorthEast')
% axis([min(sigma_t.^2),max(sigma_t.^2),0,1])
grid on
box on
str=['Alpha_1N_p','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_Alpha_2N.mat
figure(2);
hold on
plot(sigma_t.^2,m_alpha,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_ML,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_P,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_S,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_T,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','KA-P','KA-S','KA-T');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','NorthEast')
% axis([min(sigma_t.^2),max(sigma_t.^2),0,1])
grid on
box on
str=['Alpha_2N_p','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load g_Alpha_1N.mat
figure(3);
hold on
plot(sigma_t.^2,m_alpha,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_ML,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_P,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_S,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_T,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','KA-P','KA-S','KA-T');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','NorthEast')
% axis([min(sigma_t.^2),max(sigma_t.^2),0,1])
grid on
box on
str=['Alpha_1N_g','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load g_Alpha_2N.mat
figure(4);
hold on
plot(sigma_t.^2,m_alpha,'c-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_ML,'m-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_P,'g->','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_S,'b-o','linewidth',linewide1,'MarkerSize',mkft)
plot(sigma_t.^2,m_alpha_T,'r-*','linewidth',linewide1,'MarkerSize',mkft)
h_leg = legend('CC','ML','KA-P','KA-S','KA-T');
xlabel('\sigma^2','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('\alpha','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[300 0 1200 1000])
set(h_leg,'Location','NorthEast')
% axis([min(sigma_t.^2),max(sigma_t.^2)+0.05,0,1])
grid on
box on
str=['Alpha_2N_g','.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径