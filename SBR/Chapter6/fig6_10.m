%%发射和接收方向图
clc;clear;close all
lambda = 2;
%%波束指向
El=90;%%俯仰
Az=90;%%方位
%%发射
N_az = 384;
d_az = 1.08;
%%接收
% N_az = 12;
% d_az = 1.08;
%%%
N_el = 12;
d_el = 1.39;
El0=0:0.5:180;
Az0=0:0.5:180;
A = zeros(length(El0),length(Az0));
for Eli=1:1:length(El0)
    Eli
    for Azi=1:1:length(Az0)
        A(Eli,Azi)=fun_Ax(N_az, N_el, d_az, d_el, lambda, Az, El, Az0(Azi), El0(Eli));
    end
end
A=abs(A);
max(max(A))
A=A/max(max(A));
[hang,lie]=find(A==max(max(A)));
[Azx,Elx]=meshgrid(Az0,El0);
figure()
mesh(Azx,Elx,db(A));
xlabel('方位向')
ylabel('俯仰向')
figure()
plot(El0,db(A(:,lie)));
xlabel('俯仰向')
grid on
figure()
plot(Az0,db(A(hang,:)));
xlabel('方位向')
grid on