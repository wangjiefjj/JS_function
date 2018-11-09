clc
clear
close all
load 不加权脉压后的t38数据.mat
[Np,Na]=size(data_pc);
L = length(data_pc{1,1});
Tensor_Pc = zeros(Na,Np,L);%%脉冲，距离单元，阵元
for i = 1:Np
    for j = 1:Na
        Tensor_Pc(j,i,:)=data_pc{i,j};
    end
end
