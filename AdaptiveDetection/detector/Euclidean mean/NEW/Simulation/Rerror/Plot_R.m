clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 10;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:4
    load p_Rerror_1N.mat
    y1=[m_errorRSCM(i),m_errorRCC(i),m_errorRML(i),m_errorRECCP(i),...
        m_errorRECCS(i),m_errorRECCT(i)];
    load p_Rerror_2N.mat
    y2=[m_errorRSCM(i),m_errorRCC(i),m_errorRML(i),m_errorRECCP(i),...
        m_errorRECCS(i),m_errorRECCT(i)];
    y = [y1;y2];
    y=y';
    x = {'SCM','CC','ML','KA-P','KA-S','KA-T'};
    %%
    figure(i)
    bar(y,1.5)
%     caxis([0,0.6])
    set(gca, 'XTickLabel', x);
    h_leg = legend('K=1N','K=2N');
    xlabel('Estimators')
    ylabel('NFN')
    set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
    set(gcf,'Position',[300 0 1200 1000])
    grid on
    box on
    str=['Rerror_p',num2str(i),'.eps'];
    print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%     saveas(gcf,str,'psc2')
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:4
    load g_Rerror_1N.mat
    y1=[m_errorRSCM(i),m_errorRCC(i),m_errorRML(i),m_errorRECCP(i),...
        m_errorRECCS(i),m_errorRECCT(i)];
    load g_Rerror_2N.mat
    y2=[m_errorRSCM(i),m_errorRCC(i),m_errorRML(i),m_errorRECCP(i),...
        m_errorRECCS(i),m_errorRECCT(i)];
    y = [y1;y2];
    y=y';
    x = {'SCM','CC','ML','KA-P','KA-S','KA-T'};
    %%
    figure(i+4)
    bar(y,1.5)
    set(gca, 'XTickLabel', x);
    set(gca, 'YLim',[0,0.6])
    h_leg = legend('K=1N','K=2N');
    xlabel('Estimators')
    ylabel('NFN')
    set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
    set(gcf,'Position',[300 0 1200 1000])
    grid on
    box on
    str=['Rerror_g',num2str(i),'.eps'];
    print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%     saveas(gcf,str,'psc2')
end