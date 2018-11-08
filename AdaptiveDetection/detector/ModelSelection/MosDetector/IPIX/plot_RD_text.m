clc
clear
close all
FontSize = 20;
markersize =10;
linewidth = 2;
%%%%%%%%%%%%IPIXæ‡¿Î∂‡∆’¿’Õº%%%%%%%%%%%%
load ('19980204_224024_IPIX.mat')%19980204_224024_IPIX%19980223_170435_IPIX
figure()
mesh(abs(sig))
[M,L] = size(sig);
if M<L
    sig = sig.';
    [M,L] = size(sig);
end
LL = 1:L;
MM = linspace(-0.5,0.5,M);
[X,Y]=meshgrid(LL,MM);
MTD = abs(fftshift(fft(sig,[],1)));
[x,y] = max(MTD);
[~,max_Range] = max(x);
figure()
plot3(max_Range,MM(y(max_Range)),1,'ro','markersize',10,'MarkerFaceColor','r')
str = [' Range cell: %.f \n Normalized Doppler: %.6f \n Normalized amplitude: %.3f'];
text(max_Range,MM(y(max_Range)),1,sprintf(str,max_Range,MM(y(max_Range)),1),...
    'VerticalAlignment','bottom','FontSize',20)
hold on
mesh(X,Y,MTD/max(max(MTD)));
ylabel('Normalized  Doppler','FontSize',FontSize)
xlabel('Range cell','FontSize',FontSize)
zlabel('Normalized amplitude','FontSize',FontSize)
set(gca,'FontSize',FontSize)
set(gcf,'Position',[700 0 1200 1000])
grid on
box on
% figure
% for i  =1:100
%     [a,t,f] = tfrstft(sig(i,:).');
%     % f=f+0.5;
%     % f=circshift(f,18);
%     [Y,X]=meshgrid(t,f);
%     subplot(2,1,1)
%     mesh(abs(a))
%     title(num2str(i))
%     subplot(2,1,2)
%     plot(abs(sig(i,:).'))
%     pause(0.5)
%     
% end
% axis([-0.5,0.5,1,28])
% figure
% plot(abs(MTD(10001,:)))
