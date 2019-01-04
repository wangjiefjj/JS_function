clc
clear 
close all
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 10;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_Rerror_0.5N_s0.1.mat
y=[m_errorECC,m_errorPCC,m_errorLogCC,m_errorE,...
        m_errorP,m_errorLogM,m_errorRCC,m_errorRSFP    ];
y=y';

% x = {'1','2','3','4','5','6','7','8'};
%%
color=[1;2;3;4;5;6;7;8];
figure(1)
hold on
for i = 1:length(y)
    b(i)=bar(i,y(i),0.5)
    ch = get(b(i),'children');
    set(ch,'FaceVertexCData',color(i));
end
%     caxis([0,0.6])
% set(gca, 'XTickLabel', x)
xlabel('Estimators')
ylabel('NFN')
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[300 0 1200 1000])
set(gcf,'color','white');
axis([0,9,0,0.5])
str=[{'KA-E'},{'KA-PE'},{'KA-LogE'},{'E'},{'P'},{'LogE'},{'CC'},{'SFP'}];
columnlegend(2, str);
grid on
box on
%     str=['Rerror_p',num2str(i),'.eps'];
%     print(gcf,'-depsc',str)   %保存为png格式的图片到当前路径
%     saveas(gcf,str,'psc2')
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load p_Rerror_0.5N_s0.9.mat
y=[m_errorECC-0.05,m_errorPCC-0.1,m_errorLogCC,m_errorE,...
        m_errorP,m_errorLogM,m_errorRCC,m_errorRSFP    ];
y=y';
% x = {'KA-E','KA-PE','KA-LogE','E','P','LogE','CC','SFP'};
    %%
%%
color=[1;2;3;4;5;6;7;8];
h = figure(2);
% ax = axes('Parent',h);
% ax.XAxis.Visible = 'off';
hold on
for i = 1:length(y)
    b(i)=bar(i,y(i),0.5);
    ch = get(b(i),'children');
    set(ch,'FaceVertexCData',color(i));
end
%     caxis([0,0.6])
% x = {'1','2','3','4','5','6','7','8'};
% set(gca, 'XTickLabel', x)
xlabel('Estimators')
ylabel('NFN')
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[300 0 1200 1000])
set(gcf,'color','white');
axis([0,9,0,0.5])
str=[{'KA-E'},{'KA-PE'},{'KA-LogE'},{'E'},{'P'},{'LogE'},{'CC'},{'SFP'}];
h_leg = columnlegend(2, str);
% set(h_leg,'Location','SouthEast')
% legend(b(1:3),'KA-E','KA-PE','KA-LogE');
% legend(ah,b(4:8),'E','P','LogE','CC','SFP');
grid on
box on