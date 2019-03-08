clc
clear 
close all
Range = 500e3;
Ns = 16;
Np = 16;
d_az = 1.08;
fps = 0;%0;
fd = -0.5:0.01:0.5;
%%
R00 =fun_SBRR( 0, 0,Range);%fun_SBRR
PSD00 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD00(i,j) = abs(s'/R00*s);
    end
end
PSD00 = abs(PSD00)/max(max(abs(PSD00)));

if length(fps)==1
    figure()
    plot(fd,db(PSD00),'r--');
    hold on 
else
    figure()
    [X,Y]=meshgrid(fd,fps);
    mesh(X,Y,db(PSD00))
    xlabel('fd')
    ylabel('fps')
end
%%
R10 =fun_SBRR( 1, 0,Range);%fun_SBRR
PSD10 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD10(i,j) = s'/R10*s;
    end
end
PSD10 = abs(PSD10)/max(max(abs(PSD00)));

if length(fps)==1
    plot(fd,db(PSD10),'k');
else
    figure()
    [X,Y]=meshgrid(fd,fps);
    mesh(X,Y,db(abs(PSD10)))
    xlabel('fd')
    ylabel('fps')
end  
%%
R01 =fun_SBRR( 0, 1,Range);%fun_SBRR
PSD01 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD01(i,j) = s'/R01*s;
    end
end
PSD01 = abs(PSD01)/max(max(abs(PSD00)));
if length(fps)==1
    plot(fd,db(PSD01),'g--');
else
    figure()
    [X,Y]=meshgrid(fd,fps);
    mesh(X,Y,db(PSD01))
    xlabel('fd')
    ylabel('fps')
end  
%%
R11 =fun_SBRR( 1, 1,Range);%fun_SBRR
PSD11 = zeros(length(fps),length(fd));
for i = 1:length(fps)
    as = exp(1j*pi*d_az*(0:Ns-1)*fps(i))/sqrt(Ns);
    for j = 1:length(fd)
        ad = exp(1j*pi*(0:Np-1)*fd(j))/sqrt(Np);
        s = kron(ad,as).';
        PSD11(i,j) = s'/R11*s;
    end
end
PSD11 = abs(PSD11)/max(max(abs(PSD00)));

if length(fps)==1
    plot(fd,db(PSD11),'b');%,max(max(abs(PSD00)))
    xlabel('fd')
    legend('无混叠无自转','有混叠无自转','无混叠有自转','有混叠有自转')
else
    figure()
    [X,Y]=meshgrid(fd,fps);
    mesh(X,Y,db(PSD11))
    xlabel('fd')
    ylabel('fps')
end  