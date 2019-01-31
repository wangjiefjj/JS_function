clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Llist = 0.3:0.3:0.99;
load 3D_1K.mat
[X,Y] = meshgrid(SNRout,1-cos2);
figure(1);
hold on

[C,h]=contour(X,Y,Pd_PGLRT_mc','k-','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)    
[C,h]=contour(X,Y,Pd_PRAO_mc','r--','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)      
[C,h]=contour(X,Y,Pd_PWALD_mc','b-.','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)      
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','NorthWest')
axis([20,40,0,1])
grid on
box on
str=['Mis3D_1K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 3D_2K.mat
[X,Y] = meshgrid(SNRout,1-cos2);
figure(2);
hold on
[C,h]=contour(X,Y,Pd_PGLRT_mc','k-','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)    
[C,h]=contour(X,Y,Pd_PRAO_mc','r--','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)      
[C,h]=contour(X,Y,Pd_PWALD_mc','b-.','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','NorthWest')
axis([20,40,0,1])
grid on
box on
str=['Mis3D_2K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 3D_4K.mat
[X,Y] = meshgrid(SNRout,1-cos2);
figure(3);
hold on
[C,h]=contour(X,Y,Pd_PGLRT_mc','k-','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)    
[C,h]=contour(X,Y,Pd_PRAO_mc','r--','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)      
[C,h]=contour(X,Y,Pd_PWALD_mc','b-.','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','NorthWest')
axis([20,40,0,1])
grid on
box on
str=['Mis3D_4K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 3D_10K.mat
[X,Y] = meshgrid(SNRout,1-cos2);
figure(4);
hold on
[C,h]=contour(X,Y,Pd_PGLRT_mc','k-','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)    
[C,h]=contour(X,Y,Pd_PRAO_mc','r--','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)      
[C,h]=contour(X,Y,Pd_PWALD_mc','b-.','linewidth',linewide1,'ShowText','on',...
        'LevelList',Llist);
clabel(C,h,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
h_leg=legend('P-SGLRT','P-SRAO','P-SWALD');
xlabel('SNR/dB','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
ylabel('cos^2{\phi}','FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gca,'FontSize',labeltsize)
set(gcf,'Position',[500 0 1200 1000])
set(h_leg,'Location','SouthWest')
axis([20,40,0,1])
grid on
box on
str=['Mis3D_10K.eps'];
print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径