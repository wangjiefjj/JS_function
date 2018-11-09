%% ������ѹ�������     �����������    �Ӳ���ʱ��
clear all;close all;clc;
% load ��������ѹ���t38����.mat
load ����Ȩ��ѹ���t38����.mat
[p_num,c_num] = size(data_pc);       % ������   ��Ԫ��
range_num = length(data_pc{1,1});   % ��������
fd_num = 2^ceil(log2(p_num));         % ������ͨ����
%% ����--������ͼ
rd_spe = zeros(range_num,fd_num); % ���ڼ�¼RD��
for cc = 1:c_num
    temp_data = [];
    for pp = 1:p_num
        temp_data = [temp_data data_pc{pp,cc}];
    end
    %%ÿ��ͨ����rd���������
    rd_spe = rd_spe+fftshift(fft(temp_data,fd_num,2),2);
end
rd_spe_dB = db(rd_spe)./2;
range = 1:range_num;
fd = linspace(-0.5,0.5,fd_num);
figure,mesh(fd,range,rd_spe_dB);xlabel('������ͨ��'),ylabel('������');title('��ͨ����ͺ�ľ���-���������');view([0 90]);colorbar;

%% �Ӳ���ʱ�׷���
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
        st = exp(2j*pi*(0:p_num-1).'*fd(mm)); % ͨ���������ն���Ϊ        fd = -2V/lambda  ����ǰ����˸���
        sss = kron(st,ss);
        Psp(nn,mm) = 1/abs(sss'*Rin*sss);
    end
end
Psp = Psp./max(Psp(:));
Psp_dB = db(Psp)./2;
figure,mesh(fd,fs,Psp_dB);xlabel('fd'),ylabel('fs');view([0 90]);colorbar,title('�Ӳ���');
figure,surf(fd,fs,Psp_dB);xlabel('fd'),ylabel('fs');view([0 90]);colorbar;title('�Ӳ���');