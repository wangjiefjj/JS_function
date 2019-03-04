clc;clear;close all
lambda = 2;
El0=90;
Az0=90;
N_az = 384;
d_az = 1.08;
% N_az = 32;
% d_az = 13.4;
N_el = 12;
d_el = 1.39;
El=0:1:180;
Az=0:0.1:180;
A = zeros(length(El),length(Az));
for Eli=1:1:length(El)
    Eli
    for Azi=1:1:length(Az)
        A(Eli,Azi)=fun_Ax( N_az, N_el, d_az, d_el, lambda, Az0, El0, Az(Azi), El(Eli));
    end
end
A=abs(A)/max(max(A));
[hang,lie]=find(A==max(max(A)));
[Azx,Elx]=meshgrid(Az,El);
figure()
mesh(Azx,Elx,db(A));
xlabel('方位向')
ylabel('俯仰向')
figure()
plot(El,db(A(:,lie)));
xlabel('俯仰向')
figure()
plot(Az,db(A(hang,:)));
xlabel('方位向')