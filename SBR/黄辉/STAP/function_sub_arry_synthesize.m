function T_synthesize= function_sub_arry_synthesize(N,Ns,M,Ms,K,Ss_pitch_sig,Ss_azimuth_sig,channel_error)

%Ns:��λά�ϳ��������
%Ms:����ά�ϳ��������

%% �кϳ�
Ts_e = [];
a = [Ss_pitch_sig(1:M-Ms+1);zeros(Ms-1,1)];
for m = 1:Ms
    Ts_e(:,m) = circshift(a,[m-1 0]);%circshift(A,[a,b]),��ʾ����a�У�����b�У�����ά����ϳɱ任����
end
Ts_a = eye(N);
T_col = kron(Ts_a,Ts_e);%��������ϳɱ任����
    % Ts = kron(Ts_a,Ts_e);%��������ϳɱ任����
    % Ts = Ts*diag(channel_vectr);
    % Tt = eye(K);
    % T_col = kron(Tt,Ts);
    
%% �кϳ�
Ts_e = [];
a = [Ss_azimuth_sig(1:N-Ns+1);zeros(Ns-1,1)];
for m = 1:Ns
    Ts_e(:,m) = circshift(a,[m-1 0]);%circshift(A,[a,b]),��ʾ����a�У�����b�У�����ά����ϳɱ任����
end
Ts_a = eye(Ms);
T_row = kron(Ts_e,Ts_a);%��������ϳɱ任����
    % Ts = kron(Ts_e,Ts_a);%��������ϳɱ任����
    % Ts = Ts*diag(channel_vectr);
    % Tt = eye(K);
    %T_row = kron(Tt,Ts);

%% �ܵ�����ϳɾ���
T_synthesize_0=T_row*T_col ;

%���ͨ�����
channel_vectr=(1+randn(Ns*Ms,1)*channel_error).*exp(j*pi*randn(Ns*Ms,1)*channel_error);%ͨ�����
T_synthesize_0=T_synthesize_0*diag(channel_vectr);
Tt = eye(K);

T_synthesize=kron(Tt,T_synthesize_0);