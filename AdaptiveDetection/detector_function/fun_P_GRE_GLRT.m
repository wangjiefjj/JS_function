function [ Tpglrt ] = fun_P_GRE_GLRT( Train,x0,H )
%FUN_PGRE 此处显示有关此函数的摘要
%   此处显示详细说明
%%GRE关系下的Persymmetric  GLRT 检测器 
[N,L] = size(Train);
%%%产生置换矩阵，反对角线都为1其余为0
J = zeros(N,N);
for i = 1:N
    J(i,N-i+1) = 1;
end
S = Train*Train';
Sp = 0.5*(S + J*conj(S)*J);
iSp = inv(Sp);
Y = [0.5*(x0+J*conj(x0)),0.5*(x0-J*conj(x0))];
ap = (Y'*iSp*H)/(H'*iSp*H);
%%<A Persymmetric GLRT for Adaptive Detection in Partially-Homogeneous Environment>
% Fai0 = Y'*iSp*Y;
% Fai1 = Y'*iSp*Y - (Y'*iSp*p*p'*iSp*Y)/(p'*iSp*p);
% s0 = det(eye(size(Fai0))+Fai0);
% s1 = det(eye(size(Fai1))+Fai1);
% %%<Shang 2018 Multichannel>%%%%%%%%%%%%%%%%%%%%%
% s1 = 1+x0'*iSp*x0-x0'*iSp*p*inv(p'*iSp*p)*p'*iSp*x0;
% s0 = 1+(x0'*iSp*x0);
%%<GRE关系下的Persymmetric  GLRT 检测器 >%%%%%%%%%%%%%%%%%%%%%
%%1
s0 = det(eye(2)+Y'*iSp*Y);
s1 = det(eye(2)+Y'*iSp*Y-Y'*iSp*H*inv(H'*iSp*H)*H'*iSp*Y);
%%2
% s0 = det(iSp*S+x0'*iSp*x0);
% s1 = det(iSp*S+x0'*iSp*x0-x0'*iSp*p*inv(p'*iSp*p)*p'*iSp*x0);
%%%%%%%%
% Pp = Rp'*iSp*p*inv(p'*iSp*p)*p'*iSp*Rp;
% s0 = det(eye(2)+Rp'*iSp*Rp);
% s1 = det(eye(2)+(Rp'*iSp*Rp-Pp));
Tpglrt = abs(s0)/abs(s1);
end

