function [MTI_FFT] = fun_MTI(Pc)
%% MTI 
h = [1 -2 1]';  % 3脉冲MTI滤波器的H(z) = 1 - 2z.^(-1) + z.^2
[M,L]=size(Pc);
for i=1:L
    Pcm(:,i) = conv(h,Pc(:,i)); %对于相同的距离单元，在脉冲间做运算：Vout = V(i) - 2V(i -1) + V(i - 2) 是这样干的
end %
[MyPcm,NyPcm]=size(Pcm);
% Lfft = M;
% Lfft = 2^(nextpow2(M));
Lfft = 2^(nextpow2(M)+3);
YPM=fft(Pcm.*(hamming(MyPcm)*ones(1,NyPcm)),Lfft); 
MTI_FFT= fftshift(YPM.*conj(YPM),1);
end

