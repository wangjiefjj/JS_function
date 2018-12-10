%% ���źŵ���ʱ����غ���
% ���룺     x��������ʽ���ź�
% �����     N��Nά������غ�����NΪ�ź�x�ĳ���
function RXC_out = RXC(x)

if (nargin == 0)
    error('At least one parameter required');
end;

[N,xcol] = size(x);


if (xcol==0)||(xcol>2),
    error('X must be a column vector');
end;

RXC_out= zeros (N,N);
t=1:N;
for n=1:N,
    ti= t(icol); taumax=min([ti-1,N-ti,round(N/2)-1]);
    tau=-taumax:taumax; indices= rem(N/2+tau,N)+1;
    RXC_out(indices,n) = x(ti+tau) .* exp(-1j*angle(x(ti-tau)));
end;

end