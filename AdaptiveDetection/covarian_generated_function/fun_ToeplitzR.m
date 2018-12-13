function [R] = fun_ToeplitzR(N,fc,sigmaf)
%FUN_T �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%     fc   �Ӳ�����Ƶ��
%     sigmaf  %%�Ӳ���չ��
%     N ά��
nn = (0:N-1)';
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
R = toeplitz(rc);
end

