function [ Point ] = fun_findPoint( cmj )
%FUN_FINDPOINT 此处显示有关此函数的摘要
%   此处显示详细说明
%找到锥角相等的点（EL，AZ）
L = 180;
EL = linspace(0,89.5,L);
AZ = linspace(0,179,L);
count = 0;
for i = 1:L
    for j = 1:L
        if abs(cmj - sin(EL(i)/180*pi) * cos(AZ(j)/180*pi))<1e-5
            count = count + 1;
            Point(count,:) =[EL(i),AZ(j)]; 
        end
    end
end
end

