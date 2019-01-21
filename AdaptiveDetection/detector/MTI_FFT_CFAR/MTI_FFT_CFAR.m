clc
clear 
close all
% 
%%
Na=8;     
Np=8;     
N=Na*Np;
optc = 'g';
opt_train = 1;%%1:SIRP,2:部分均匀
L=round(2*N); 
cos2=1;%%%失配
PFA=1e-4;% PFA=1e-4;
%%各种比
SNRout = 20; % 输出信噪比SNR
CNRout = 15; %杂噪比
JNRout = 15; %干噪比
SNRnum=10.^(SNRout/10);
CNRnum=10.^(CNRout/10);
JNRnum=10.^(JNRout/10);
%%虚警率和MC次数
MonteCarloPfa=round(1/PFA*100);
MonteCarloPd=1e4;
%%信号
fs=[0.2];
theta = ones(length(fs),1);
nn=(0:N-1)';
vt = exp(1i*2*pi*nn*fs)/sqrt(N); %%%%%% 导向矢量
% [UU,SS,VV]=svd(iR*H);
% vt_v=UU(:,end); %导向矢量和杂波的正交分量vt^H*iR*vt_v==0,
% q = 2*vt_v;%%GRE中的q
% R0 = R + q*q';
%% 各种协方差
%%杂波协方差
fc = 0;
sigmaf = 0.03; %%杂波谱展宽
% rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
% Rc = CNRnum * toeplitz(rc);
Rc = CNRnum * fun_rho(0.9,N,1,0.0);  
%%干扰协方差
fj = -0.2;
jam = exp(-1i*2*pi*nn*fj);
Rj = JNRnum * jam*jam';
%%检测单元协方差test
Rt = eye(N) + Rc + Rj;
%%参考单元协方差secondary 
Rs = eye(N) + Rc;
%% 回波脉压的数据
Train = fun_TrainData(optc,N,L,Rs,3,1,opt_train);
x = fun_TrainData(optc,N,1,Rs ,3,1,opt_train);
alpha=sqrt(SNRnum/abs(vt'/Rs*vt)); 
x=alpha*vt+x;%+pp;    
Pc = [Train(:,1:L/2),x,Train(:,L/2+1:end)];%%包含目标的脉压结果
RD=fftshift(fft(Pc,[],1),1);
fd = linspace(-0.5,0.5,N);
[X,Y]=meshgrid(1:L+1,fd);
mesh(X,Y,db(abs(RD),'power'))
%% MTI 
h = [1 -2 1]';  % 3脉冲MTI滤波器的H(z) = 1 - 2z.^(-1) + z.^2
for i=1:L+1
    Pcm(:,i) = conv(h,Pc(:,i)); %对于相同的距离单元，在脉冲间做运算：Vout = V(i) - 2V(i -1) + V(i - 2) 是这样干的
end %
[MyPcm,NyPcm]=size(Pcm);
% Lfft = N;
% Lfft = 2^(nextpow2(N));
Lfft = 2^(nextpow2(N)+3);
YPM=fft(Pcm.*(hamming(MyPcm)*ones(1,NyPcm)),Lfft); 
YPM= fftshift(YPM.*conj(YPM),1);
YPMdB=db(abs(YPM),'power');
fd = linspace(-0.5,0.5,Lfft);
[X,Y]=meshgrid(1:L+1,fd);
figure()
mesh(X,Y,YPMdB)
levels=(max(YPMdB(:))+[-1 -5 -10 -15 -20 -25]);
figure
contour(fd,1:L+1,YPMdB',levels)
%% CFAR,18个单元的CA-CFAR
Reference = 18;
cfar = ones(Reference+5,1)/Reference;
CUT = Reference/2+3;
cfar(CUT-2:CUT+2)=0; % Cell Under Test单元的index是12，它左、右各有 2个guard cells，9个averaging cells
Mean = zeros(size(YPM)); %%各检测单元的参考单元均值
for i=1:Lfft
    temp = conv(YPM(i,:),cfar); 
    Mean(i,:) = temp(12:end-11);
end
MeandB=db(abs(Mean),'power');
alpha=(Reference)*(PFA^(-1/(Reference))-1);%求alpha
threshold=alpha*Mean;%%求门限;
detected = YPM.*(YPM>threshold);
figure
mesh(X,Y,abs(detected))