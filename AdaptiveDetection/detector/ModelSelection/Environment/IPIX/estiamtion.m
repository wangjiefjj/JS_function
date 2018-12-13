clc
clear 
close all
%%%%%%¹À¼Æ²ÎÊý
% Read_Display_Data
Data_process
load(matFile)
tic
t = [9:29];
offset = 50000;%50000;
parfor i = 1:10000%nsweep
    r = sig(i+offset,t);
    [lambda(i),mu(i)] = fun_IG_ML(r);
end
lambda_m=mean(lambda)
mu_m=mean(mu)
toc

% MTD = abs(fftshift(fft(sig,[],1)));
% mesh(MTD)
% save(matFile,'R','Range','M','N','sig','lambda','mu');