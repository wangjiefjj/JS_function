clc
clear 
close all
% 根据文献《A New CFAR Detection Test for Radar》编写的
% 自�?应匹配滤波（AMF�?
%%暂时认为每个脉冲的到达角相同�?
%%%增益是MN*脉压增益，M=1，N=1,信号幅度�?�?
%%%SNR信噪比对应的�?��概率=SNR-10*log10（MN*脉压增益）时的检测概�?
fc = 1e9;       %载频
C = 3e8;        %光�?
lamda = C/fc;   %波长
tao = 127e-6;    %脉宽
B = 1e6;        %带宽
mu = B/tao;     %调频�?
Fs = 2*B;       %采样频率
Ts = 1/Fs;
t = -tao/2:Ts:tao/2-Ts;  %快时�?
L = length(t);
R = 0e3;
dt = 2*R/C;
M = 2;                  %脉冲�?
%%阵列参数
N = 1;                  %阵元�?
d = 0.5*lamda;          %阵元间隔
theta = 0;%jiao2hu(10);    %到达�?
St = exp(1j*2*pi*d*(0:N-1).'/lamda*sin(theta));%导向矢量
S = St;
S = repmat(S,[M,1]);    %导向矢量
%%%接收处理(脉压)
signal = exp(-1j*2*pi*(fc*t+0.5*mu*t.^2) );
h_fft = fft(signal);
MC = 1000;
%%%
Pfa = 1e-6;     %虚警�?
r0 = -log(Pfa); %门限
SNR = -17:0.1:0;
for SNR_i = 1:length(SNR)
    SNR_i
    num = 0;
        for MC_i = 1:MC
            Pc = zeros(N,M*L);
%             A = 1*rand()+1j*1*rand();
            A=1;
            for i =1:M %脉冲
                for j =1:N %阵列
                    echo = A*exp(1j*2*pi*(fc*dt+0.5*mu*(t-dt).^2) );
                    echo = awgn(echo,SNR(SNR_i),'measured');
                    echo_fft = fft(echo);
                    Pct = ifftshift(ifft(h_fft.*echo_fft)).*repmat(St(j),[1,L]);
                    Pct = Pct.*repmat(St(j),[1,L]);
                    Pc(j,1+(i-1)*L:i*L) = Pct;
                end
            end
          detect_R=L/2+1;
    %     detect_R=round(L/4);%%没目标的位置
            x1 = [];
            x2 = [];
            x = [];
            y = [];
            for i =1:M
                x1 = Pc(:,1+L*(i-1):detect_R+L*(i-1)-10);
                x2=  Pc(:,detect_R+L*(i-1)+10:L*i);
                x = cat(1,x,cat(2,x1,x2));
                y = cat(1,y,Pc(:,detect_R+L*(i-1)));
            end
            M_gu = x*x'/length(x);
            r = abs((abs(S.'*inv(M_gu)*y))^2/(S.'*inv(M_gu)*S));%GLRT—AMF�?���?
            if r > r0
                num = num+1;
            end
        end
        proportion(SNR_i) = num/MC;
end
plot(SNR,proportion)

% figure()
% plot(abs(Pc(1,:)))
% %%%理论�?
% %%%
Pfa = 1e-6;     %虚警�?
r0 = -log(Pfa); %门限
SNR = 0:0.1:20;
alpha = SNR2real(SNR,10);
n = 0:150;
for i =1:length(alpha)
    Pd(i) = sum(alpha(i).^n.*exp(-alpha(i))./factorial(n).*gammainc(n+1,r0));
end
% figure()
hold on
plot(SNR,Pd,'r')


