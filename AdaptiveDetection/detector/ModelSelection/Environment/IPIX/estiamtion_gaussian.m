clc
clear 
close all
%%%%%%���Ʋ���
% Read_Display_Data
tic
N=10000;
L=21;
R = fun_rho(0.9,N,1,0.1);
[Train,tauk] = fun_TrainData('g',N,L,R);%%������ѵ������,Э�������ΪrouR�ĸ�˹�Ӳ�
parfor i = 1:10000%nsweep
    [lambda(i),mu(i)] = fun_IG_ML(Train(i,:));
end
mean(lambda)
mean(mu)
toc

% MTD = abs(fftshift(fft(sig,[],1)));
% mesh(MTD)
% save(matFile,'R','Range','M','N','sig','lambda','mu');