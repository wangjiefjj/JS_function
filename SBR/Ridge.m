clc
clear 
close all
Range = 500e3;
Ns = 16;
Np = 16;
d_az = 1;
fps = -0.5:0.01:0.5;%0;
fd = -0.5:0.01:0.5;
%%
R00 =fun_SBRR( 0, 0,Range);%fun_SBRR
PSD00 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    i
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD00(i,j) = abs(s'*R00*s)^2;
    end
end
PSD00 = abs(PSD00)/max(max(abs(PSD00)));
figure()
[X,Y]=meshgrid(fd,fps);
imagesc(fd,fps,db(PSD00))
%     mesh(X,Y,db(PSD00))
xlabel('fd')
ylabel('fps')
%%
R10 =fun_SBRR( 1, 0,Range);%fun_SBRR
PSD10 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD10(i,j) = abs(s'*R10*s)^2;
    end
end
PSD10 = abs(PSD10)/max(max(abs(PSD10)));

figure()
[X,Y]=meshgrid(fd,fps);
imagesc(fd,fps,db(PSD10))
%     mesh(X,Y,db(abs(PSD10)))
xlabel('fd')
ylabel('fps') 
%%
R01 =fun_SBRR( 0, 1,Range);%fun_SBRR
PSD01 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD01(i,j) = abs(s'*R01*s)^2;
    end
end
PSD01 = abs(PSD01)/max(max(abs(PSD01)));
figure()
[X,Y]=meshgrid(fd,fps);
imagesc(fd,fps,db(PSD01))
%     mesh(X,Y,db(PSD01))
xlabel('fd')
ylabel('fps')  
%%
R11 =fun_SBRR( 1, 1,Range);%fun_SBRR
PSD11 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD11(i,j) = abs(s'*R11*s)^2;
    end
end
PSD11 = abs(PSD11)/max(max(abs(PSD11)));
figure()
[X,Y]=meshgrid(fd,fps);
imagesc(fd,fps,db(PSD11))
%     mesh(X,Y,db(PSD11))
xlabel('fd')
ylabel('fps') 