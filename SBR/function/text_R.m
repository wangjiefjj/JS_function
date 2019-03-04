clc
clear 
close all
count = 0;
for range = 550e3:10e3:1000e3
    count = count+1
    R1 =fun_SBRR( 1, 1,range);
    rankR1(count)=rank(R1);
    R0 =fun_SBRR( 1, 0,range);
    rankR0(count)=rank(R0);
end
% R1 =fun_SBRR( 1, 1,550e3);
% E=abs(eig(R1));
% E = db(sort(E,'descend')).';
% plot(E,'b','LineWidth',2)
% hold on
% R1 =fun_SBRR( 1, 0,550e3);
% E=abs(eig(R1));
% E = db(sort(E,'descend')).';
% plot(E,'r','LineWidth',2)