clc
clear
close all
load ����Ȩ��ѹ���t38����.mat
[Np,Na]=size(data_pc);
L = length(data_pc{1,1});
Tensor_Pc = zeros(Na,Np,L);%%���壬���뵥Ԫ����Ԫ
for i = 1:Np
    for j = 1:Na
        Tensor_Pc(j,i,:)=data_pc{i,j};
    end
end
