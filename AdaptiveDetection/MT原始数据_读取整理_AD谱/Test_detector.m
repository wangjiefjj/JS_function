clc
clear 
close all
load 不加权脉压后的t38数据全维数据整理.mat
[N,L] =  size(data_st);
c_num=14;
p_num=16;
N = c_num*p_num;
Th_ACE = fun_Th_ACE(403,224,1e-6);
fs = 0.1;
fd = 0.1;
ss = exp(2j*pi*(0:c_num-1).'*fs);
st = exp(2j*pi*(0:p_num-1).'*fd);
s = kron(st,ss);
s = s/(sqrt(N));
SNRout = -5:25;
SNRnum = 10.^(SNRout/10);
Re = data_st*data_st'./403;
iRe = inv(Re);
alpha=sqrt(SNRnum/abs(s'*iRe*s));

%%ANMF不同信噪比下的空时图
% L_st = 20;
% fdd = linspace(-0.5,0.5,L_st);
% fss = linspace(-0.5,0.5,L_st);
% for i =1:length(SNRnum)
%     i
%     x0 = data_st(:,200);
%     x0 = x0 + alpha(i)*s;
%     for nn = 1:L_st
%         ss = exp(2j*pi*(0:c_num-1).'*fss(nn));
%         for mm = 1:L_st
%             st = exp(2j*pi*(0:p_num-1).'*fdd(mm)); % 通常将多普勒定义为        fd = -2V/lambda  所以前面加了负号
%             sss = kron(st,ss)/sqrt(N);
%             ACE_st(nn,mm,i) = fun_ANMF(Re,x0,sss);
%         end
%     end
% end
% for i = 1:31
%     pause(0.5)
%     mesh(abs(ACE_st(:,:,i)))
%     str=[num2str(SNRout(i)),'dB'];
%     title(str)
% end

%%ANMF不同信噪比下的个距离单元滤波值
temp1 = data_st(:,120);
temp2 = data_st(:,110);
% R_LogM = fun_RLogEMean(data_st,4);
R_NSCM =fun_NSCMN(data_st);
for i =1:length(SNRnum)
    i
    data_st(:,120) = temp1 + alpha(i)*s;
    data_st(:,110) = temp2 + alpha(i)*s;
    count = 0;
    for nn = 100:150
        count = count +1;
        ACE_range(i,count) =  fun_ANMF(Re,data_st(:,nn),s);
        LogM_rang(i,count) = fun_ANMF(R_NSCM,data_st(:,nn),s);
    end
end
for i = 1:31
    pause(0.5)
%     plot(ACE_range(i,:),'b')
    plot(LogM_rang(i,:),'k')
    hold on 
    plot(ones(1,50)*Th_ACE,'r')
    str=[num2str(SNRout(i)),'dB'];
    title(str)
end