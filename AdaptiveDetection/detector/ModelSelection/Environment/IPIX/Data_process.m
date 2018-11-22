Read_Display_Data
load(matFile) 
N = 8;
% Range = 17;
[row,col] = size(sig);
if col>row
    sig = sig.';
end
[M,L] = size(sig);

figure()
mesh(abs(sig))
sig_test_t = sig(2:N+1,:);

[X,Y]=meshgrid(1:L,linspace(-0.5,0.5,M));
MTD = abs(fftshift(fft(sig,[],1)));
[x,y] = max(MTD);
[~,Range] = max(x);
% figure()
% mesh(X,Y,MTD);
R = zeros(N,N);
L_KA = L-1;
for i = 1:L_KA
    X_KA = sig((i-1)*N+1:i*N,Range);
    R = R+X_KA*X_KA'/L_KA;
end
% figure(3)
% mesh(abs(R_KA))
save(matFile,'R','Range','M','N','sig');
