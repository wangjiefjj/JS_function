function [ GrazeTH,R,Rs] = fun_GrazeTH(H, EL_beam,opt )
%FUN_ 此处显示有关此函数的摘要
%   此处显示详细说明
% EL_beam: 俯仰向波束宽度。rad
%% 计算远近地点掠射角
Re = 6373; %flat earth radius，km
[EL,R,Rs] = fun_ELAngle(H); 
EL = EL/180*pi;
if opt==1 %%远地点
    temp = (1+H/Re)*sin(EL+EL_beam/2);
    GrazeTH = abs(acos(temp));
else
	temp = (1+H/Re)*sin(EL-EL_beam/2);
	GrazeTH = abs(acos(temp));  
end


end

