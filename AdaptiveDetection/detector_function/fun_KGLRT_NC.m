function [ kglrt_nc ] = fun_KGLRT_NC( Train,x0,p,N)
%FUN 此处显示有关此函数的摘要
%   此处显示详细说明
%% 分组积累KGLRT
%Train： 全部训练数据
% x0：全部待检测单元
% p: 全部的导向矢量
% N: 分的组数
[M,L]=size(Train);
%% 开始
kglrt_nc = 1;
m = M/N; %每组的维数
for i = 1:N
    Train_t = Train((i-1)*m+1:i*m,:);
    S = fun_SCMN(Train_t);
    x0_t = x0((i-1)*m+1:i*m);
%     p_t = p((i-1)*m+1:i*m);
    p_t = p(1:m);
    kglrt_nc = kglrt_nc * fun_KGLRT(S,x0_t,p_t);
end
end

