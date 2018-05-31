function [ EL, R, Rs ] = fun_ELAngle( H,R)
%FUN_ELANGLE 此处显示有关此函数的摘要
%   此处显示详细说明
%% 俯仰角与地距的关系，公式（4.7）
%%
Re = 6373; %flat earth radius，km
if nargin==1
    [R,Rs]=fun_RsR(H);
    EL=(asin(Re./Rs.*sin(R./Re)));
    EL = EL/pi*180;
else
    Rs = fun_R2Rs(H,R);
    EL=(asin(Re./Rs.*sin(R./Re)));
    EL = EL/pi*180;
end


end

