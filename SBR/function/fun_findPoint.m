function [ Point ] = fun_findPoint( cmj )
%FUN_FINDPOINT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%�ҵ�׶����ȵĵ㣨EL��AZ��
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

