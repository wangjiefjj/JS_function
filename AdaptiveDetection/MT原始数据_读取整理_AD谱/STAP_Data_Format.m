%% ��Mountain Topԭʼ���������STAP���ݸ�ʽ
clear all; close all;clc;
load t38pre01v1_cpi_6.mat
% load stap3001v1.mat
cpi_num = 1; % �������д洢������CPI������ cpi_num��ȡֵΪ1��2
cpi_data = cpi1; % ��һ��CPI����
%% ��ȡ�ͼ������
c = 3e8; % ��Ų������ٶ�
fc = fxmit(cpi_num); % ��ȡ��CPI���ݷ���Ƶ��
lambda = c/fc; % ���㲨��
azimuth_point = azxmit(cpi_num); % ����ͽ��շ�λ�� ��Ӧ����������   ��λ��
delta_d = 0.3448; % ��Ԫ���
pulse_num = npulses(cpi_num); % ������
tr = pri(cpi_num); % �����ظ�����
fr = 1/tr; % �����ظ�Ƶ��
tau = tpulse(cpi_num); % ����������
% t_fs_begin = trecord(cpi_num)*1e-6; % ��ʼ�������ݵ�ʱ�䣨�ӷ������忪ʼ����,��λ ΢�룬���Գ��� 1e-6�� ����������Ӧ
t_fs_begin = 881e-6; % ������˵��881΢�뿪ʼ
ts = 1e-6; % ����������
fs = 1/ts; % �������Ƶ��
delta_range_gate = c./(2*fs);
[r1,c1] = size(cpi_data); % r1��Ӧ������������������ĳ˻�(16*403)      c1��Ӧ����ͨ����
ch_num = c1; % ͨ����
range_num = r1/pulse_num; % ��������
B_lfm = 500e3; % LFM�źŴ���
Vp = 12.2*0.0254/tr; % ƽ̨�ٶ� �����ƶ�����������λ���ĵķ���(IDPCA)��Чƽ̨���˶�������PRI�˶� 12.2 inches

%% ������ѹ��������Ȩ��ѹ
%���� LFM�ź���ѹ�ο�����
Kr = B_lfm/tau; % LFM�źŵ�Ƶб��
t = -tau/2:ts:tau/2-ts;% t = 0:ts:tau-ts;
s_ref = exp(1i*pi*Kr*t.^2); % LFM�ź���ѹ�ο��ź�
nfft = length(s_ref)+range_num-1; % fft����
S_ref = fft(fliplr(conj(s_ref)),nfft); % Ƶ��ο�����
%���� ��������
count = 0;
for cc = 1:ch_num
    for pp = 1:pulse_num
        temp = cpi_data((pp-1)*range_num+1:pp*range_num,cc); % �� cc ��ͨ������ pp ����������
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
save ����Ȩ��ѹ���t38����.mat data_pc