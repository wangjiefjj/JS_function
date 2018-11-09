%% 分析脉压后的数据     距离多普勒谱    杂波空时谱
clear all;close all;clc;
% load 汉明窗脉压后的t38数据.mat
load 不加权脉压后的t38数据.mat
[p_num,c_num] = size(data_pc);       % 脉冲数   阵元数
range_num = length(data_pc{1,1});   % 距离门数
fd_num = 2^ceil(log2(p_num));         % 多普勒通道数
%% 距离--多普勒图
rd_spe = zeros(range_num,fd_num); % 用于记录RD谱
for cc = 1:c_num
    temp_data = [];
    for pp = 1:p_num
        temp_data = [temp_data data_pc{pp,cc}];
    end
    %%每个通道的rd结果加起来
    rd_spe = rd_spe+fftshift(fft(temp_data,fd_num,2),2);
end
rd_spe_dB = db(rd_spe)./2;
range = 1:range_num;
fd = linspace(-0.5,0.5,fd_num);
figure,mesh(fd,range,rd_spe_dB);xlabel('多普勒通道'),ylabel('距离门');title('各通道求和后的距离-多普勒输出');view([0 90]);colorbar;

%% 杂波空时谱分析
data_st = [];
for pn = 1:p_num
    temp_data = [];
    for cn = 1:c_num
        temp_data = [temp_data data_pc{pn,cn}];
    end
    data_st = [data_st;temp_data.'];
end
% data_st1 = data_st(:,1:350);
% Re = data_st1*data_st1'./350;
Re = data_st*data_st'./403;
[V0,D0] = eig(Re);
d1 = diag(D0);
d2 = sort(d1,'descend');
figure,plot(db(d2)./2,'*-');
Rin = inv(Re);%+1e3*eye(c_num*p_num));

fd = linspace(-0.5,0.5,165);
fs = linspace(-0.5,0.5,165);
for nn = 1:165
    ss = exp(2j*pi*(0:c_num-1).'*fs(nn));
    for mm = 1:165
        st = exp(2j*pi*(0:p_num-1).'*fd(mm)); % 通常将多普勒定义为        fd = -2V/lambda  所以前面加了负号
        sss = kron(st,ss);
        Psp(nn,mm) = 1/abs(sss'*Rin*sss);
    end
end
Psp = Psp./max(Psp(:));
Psp_dB = db(Psp)./2;
figure,mesh(fd,fs,Psp_dB);xlabel('fd'),ylabel('fs');view([0 90]);colorbar,title('杂波谱');
figure,surf(fd,fs,Psp_dB);xlabel('fd'),ylabel('fs');view([0 90]);colorbar;title('杂波谱');