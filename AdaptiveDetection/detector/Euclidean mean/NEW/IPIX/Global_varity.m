clc
clear 
close all
%%%%%%%��������ȫ�ֱ�������ʹ��
cdfFile =  '19980223_170435_ANTSTEP.CDF';
Range = 8;
cdfFile_t = cdfFile;
if length(cdfFile_t)>17
    cdfFile_t(17:27)=[];
else
    cdfFile_t(4:7)=[];
end
matFile = [cdfFile_t,'_IPIX.mat'];