%%计算LVD

clc;clear;close all
%信号%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc1 = 100;%载频
mu1 = 70;%调频斜率
fc2 = 100;%载频
mu2 = -40;%调频斜率
fc3 = -70;%载频
mu3 = 70;%调频斜率
% B = 1e6;%带宽
Fs = 512;%采样频率
Ts = 1/Fs;
t = (-Fs/2:Fs/2-1)*Ts;
tao = (-Fs/2:Fs/2-1)*Ts;
L = length(t);%%快时间频点数
L_tao = length(tao);
f1 = linspace(-Fs/2,Fs/2-1,L);
f_u = linspace(-Fs/2,Fs/2-1,L_tao);
A1 = 1;
A2 = 1;
A3 = 1;
s1 = A1 * exp(1j*2*pi*(fc1 * t + 0.5*mu1*t.^2))';
% angle_s1 = angle(s1);
s2 = A2 * exp(1j*2*pi*(fc2 * t + 0.5*mu2*t.^2))';
s3 = A3 * exp(1j*2*pi*(fc3 * t + 0.5*mu3*t.^2))';
s = s1+s2+s3;
% s = hilbert(s);
figure()
plot(real(s))
title('信号')
s_fft = (fft(s));
figure()
plot(f1,abs(fftshift(s_fft)))
title('信号频谱')
%%%尺度变换参数
a = 0.1;
q = a/Ts;
h = 1/a;
s1 = JS_RXC3(s,q); %%变尺度瞬时自相关

%%%s1=  [(tao1,t1),(tao1,t2),(tao1,t3)
%%%      (tao2,t1),(tao2,t2),(tao2,t3)
%%%      (tao3,t1),(tao3,t2),(tao3,t3)]
figure()
% mesh(abs(s1))
imagesc(abs(s1))
title('参数对称瞬时自相关')

s1_fft = fftshift(fft(s1,[],1),1);
figure()
imagesc(t,f_u,flipud(abs(s1_fft)));
% hold off
title('参数对称瞬时自相关频谱')

% figure;plot(f_u,abs(s1_fft(:,512)))
 
% N =length(s);
% R = zeros(N,N);
% for jj = 1:N
%     taumax = min([jj-1 N-jj floor(N/2)]);
%     tau = -taumax:taumax;
%     indices = floor(N/2)+tau+1;
%     R(indices,jj) = s(jj+tau,1) .* conj(s(jj-tau,1));  
% end
% 
% figure;imagesc(abs(R))

%%变尺度FT
tic
S = JS_SFT2(s1,a,h,Ts); 
toc

figure()
% mesh(abs(S))
imagesc(abs(S))
title('变尺度FT')
LVD = fftshift(fft(S,[],2),2);
[F,Mu] = meshgrid(f1,f_u);
figure()
mesh(-Mu,-F,abs(LVD))
xlabel('载频');ylabel('调频率')
figure;
imagesc(abs(LVD))
xlabel('载频');ylabel('调频率')

% 调频率方向的分辨能力
 figure;plot((-L/2+1:L/2),flipud(abs(LVD(:,56))),'Linewidth',1)
 xlabel('调频率')
 ylabel('幅值')
 
 % 载频方向的分辨能力
 figure;plot((-L/2+1:L/2),fliplr(abs(LVD(250,:))),'Linewidth',2)
 xlabel('载频')
 ylabel('幅值')
 