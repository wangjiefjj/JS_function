clc
clear 
close all
%%%%%%���Ʋ���
% Read_Display_Data
Data_process
load(matFile)
tic
t = [9:29];
parfor i = 1:10000%nsweep
    r = sig(i,t);
    [lambda(i),mu(i)] = fun_IG_ML(r);
end
mean(lambda)
mean(mu)
toc

% MTD = abs(fftshift(fft(sig,[],1)));
% mesh(MTD)
% save(matFile,'R','Range','M','N','sig','lambda','mu');