Read_Display_Data
load(matFile) 
N = 8;
labeltsize=35;
fw = 'normal'; %%是否加粗斜体之类
fn='Times New Roman';
linewide1=3;
mkft = 15;
% Range = 17;
[row,col] = size(sig);
if col>row
    sig = sig.';
end
[M,L] = size(sig);

figure()
imagesc(abs(sig)/max(max(abs(sig))))
xlabel('Range Cell')
ylabel('Pulse Number ')
hcb=colorbar('NorthOutside');
title(hcb,'Normalized Amplitude')
set(gca,'FontSize',labeltsize,'FontWeight',fw,'FontName',fn)
set(gcf,'Position',[100 0 1200 1000])
% set(h_leg,'Location','NorthWest')
grid on
box on