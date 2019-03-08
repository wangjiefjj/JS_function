clc
clear 
close all
Ns = 16;
Np = 16;
d_az = 1.08;
fps = 0;%0;
fd = -1:0.01:1;
Range = 500e3:10e3:2400e3;
as = exp(1j*pi*d_az*(0:Ns-1)*fps)/sqrt(Ns);
%%
PSD00 = zeros(length(Range),length(fd));
PSD10 = zeros(length(Range),length(fd));
PSD01 = zeros(length(Range),length(fd));
PSD11 = zeros(length(Range),length(fd));
for i = 1:length(Range)
    i
    R00 =fun_SBRR( 0, 0,Range(i));%fun_SBRR
    R10 =fun_SBRR( 1, 0,Range(i));%fun_SBRR
    R01 =fun_SBRR( 0, 1,Range(i));%fun_SBRR
    R11 =fun_SBRR( 1, 1,Range(i));%fun_SBRR
    parfor j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD00(i,j) = s'/R00*s;
        PSD10(i,j) = s'/R10*s;
        PSD01(i,j) = s'/R01*s;
        PSD11(i,j) = s'/R11*s;
   end
end
[X,Y]=meshgrid(fd,Range/1000);
%%
figure()
imagesc(fd,Range/1000,db(abs(PSD00),max(max(abs(PSD00)))))
view(0,-90)
% mesh(X,Y,db(abs(PSD00)))
xlabel('fd')
ylabel('Range')
%%
figure()
imagesc(fd,Range/1000,db(abs(PSD10),max(max(abs(PSD00)))))
view(0,-90)
% mesh(X,Y,db(abs(PSD10)))
xlabel('fd')
ylabel('Range')
%%
figure()
imagesc(fd,Range/1000,db(abs(PSD01),max(max(abs(PSD00)))))
view(0,-90)
% mesh(X,Y,db(abs(PSD01)))
xlabel('fd')
ylabel('Range')
%%
figure()
imagesc(fd,Range/1000,db(abs(PSD11),max(max(abs(PSD00)))))
view(0,-90)
% mesh(X,Y,db(abs(PSD11)))
xlabel('fd')
ylabel('Range')
% save('RD3.mat')