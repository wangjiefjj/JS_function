function [R] = fun_ToeplitzR(N,fc,sigmaf)
%FUN_T 此处显示有关此函数的摘要
%   此处显示详细说明
%     fc   杂波中心频率
%     sigmaf  %%杂波谱展宽
%     N 维数
nn = (0:N-1)';
rc =  exp(-1i*2*pi*nn*fc-2*(nn*pi*sigmaf).^2);
R = toeplitz(rc);
end

