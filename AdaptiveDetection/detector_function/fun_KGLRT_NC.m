function [ kglrt_nc ] = fun_KGLRT_NC( Train,x0,p,N)
%FUN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% �������KGLRT
%Train�� ȫ��ѵ������
% x0��ȫ������ⵥԪ
% p: ȫ���ĵ���ʸ��
% N: �ֵ�����
[M,L]=size(Train);
%% ��ʼ
kglrt_nc = 1;
m = M/N; %ÿ���ά��
for i = 1:N
    Train_t = Train((i-1)*m+1:i*m,:);
    S = fun_SCMN(Train_t);
    x0_t = x0((i-1)*m+1:i*m);
%     p_t = p((i-1)*m+1:i*m);
    p_t = p(1:m);
    kglrt_nc = kglrt_nc * fun_KGLRT(S,x0_t,p_t);
end
end

