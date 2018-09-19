%% 将Mountain Top原始数据整理成STAP数据格式
clear all; close all;clc;
load t38pre01v1_cpi_6.mat
% load stap3001v1.mat
cpi_num = 1; % 该数据中存储了两个CPI的数据 cpi_num可取值为1和2
cpi_data = cpi1; % 第一个CPI数据
%% 读取和计算参数
c = 3e8; % 电磁波传播速度
fc = fxmit(cpi_num); % 读取该CPI数据发射频率
lambda = c/fc; % 计算波长
azimuth_point = azxmit(cpi_num); % 发射和接收方位角 对应于正北方向   单位度
delta_d = 0.3448; % 阵元间距
pulse_num = npulses(cpi_num); % 脉冲数
tr = pri(cpi_num); % 脉冲重复周期
fr = 1/tr; % 脉冲重复频率
tau = tpulse(cpi_num); % 发射脉冲宽度
% t_fs_begin = trecord(cpi_num)*1e-6; % 开始接收数据的时间（从发射脉冲开始算起,单位 微秒，所以乘了 1e-6） 用于与距离对应
t_fs_begin = 881e-6; % 资料上说从881微秒开始
ts = 1e-6; % 距离采样间隔
fs = 1/ts; % 距离采样频率
delta_range_gate = c./(2*fs);
[r1,c1] = size(cpi_data); % r1对应距离采样数与脉冲数的乘积(16*403)      c1对应接收通道数
ch_num = c1; % 通道数
range_num = r1/pulse_num; % 距离门数
B_lfm = 500e3; % LFM信号带宽
Vp = 12.2*0.0254/tr; % 平台速度 采用移动发射天线相位中心的方法(IDPCA)等效平台的运动，相邻PRI运动 12.2 inches

%% 距离脉压——不加权脉压
%—— LFM信号脉压参考函数
Kr = B_lfm/tau; % LFM信号调频斜率
t = -tau/2:ts:tau/2-ts;% t = 0:ts:tau-ts;
s_ref = exp(1i*pi*Kr*t.^2); % LFM信号脉压参考信号
nfft = length(s_ref)+range_num-1; % fft点数
S_ref = fft(fliplr(conj(s_ref)),nfft); % 频域参考函数
%—— 整理数据
count = 0;
for cc = 1:ch_num
    for pp = 1:pulse_num
        temp = cpi_data((pp-1)*range_num+1:pp*range_num,cc); % 第 cc 个通道，第 pp 个脉冲数据
%         figure,plot(db(temp));
        temp_f = fft(temp,nfft);
        temp_pcf = temp_f.*S_ref.';
        temp_pc = ifft(temp_pcf);
        temp_pc_ef = temp_pc(length(s_ref):end);
%         figure,plot(db(temp_pc_ef));
        count = count+1;
        data_pc{pp,cc} = temp_pc_ef;
    end
end
save 不加权脉压后的t38数据.mat data_pc