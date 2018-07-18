function [ Range ] = fun_JWH_range( Aj,Aw,Ah,Bj,Bw,Bh )
%FUN_JWH_RANGE 此处显示有关此函数的摘要
%   此处显示详细说明
% Aj：经度，deg
% Bw：纬度，deg
% h：距地心的距离 
Aj = Aj/180*pi;
Aw = Aw/180*pi;
Bj = Bj/180*pi;
Bw = Bw/180*pi;
Ax = Ah * cos(Aw) * cos(Aj);
Ay = Ah * cos(Aw) * sin(Aj);
Az = Ah * sin(Aw);
Bx = Bh * cos(Bw) * cos(Bj);
By = Bh * cos(Bw) * sin(Bj);
Bz = Bh * sin(Bw);
Range = sqrt((Ax-Bx)^2 + (Ay-By)^2 + (Az-Bz)^2);
end

