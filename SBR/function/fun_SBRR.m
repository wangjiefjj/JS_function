function [ Rcn,El0,d_az,lambda,fr,Np,N_az ] = fun_SBRR( isRu, isRotation,R)
%%%����P156Ҳ��6.2�Ĳ�������
if nargin <3
    R = 1300e3;
end
%% �״�ϵͳ����
fo = 1200e6; %9600e6 450e6     % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4e6;                     % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 500;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
Np = 16;  %16                  % Number of Pulses per CPI.
Tp = 200e-6;                  % Pulse Width in sec.
N_az = 16;  %16               % Number of Array Antenna Elements, ����ʱ��384,32*12
N_el = 1;
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d_az = lambda/2*1;              % ���
% d_az = lambda/2*1.08;              % ���
d_el = lambda/2*1.39;              % ���
% Azimuth angle in degrees ��λ��:
% Az = linspace(-0,179,Nc);
Az = 0:1:180;
EL = 0:1:90;
Nc = length(Az);
LAz = length(Az);
f = zeros(1,LAz);
AF = ones(1,LAz);           % Array Factor pre-allocation. ����������
Sys_DOF = Np*N_az*N_el;
dR = c/2/B;                   % Radar Range Resolution in meters.
%% Platform Parameters:
H = 506e3;%9e3 %700e3                    % Platform altitude in meters.
Rs = fun_R2Rs(H,R);
if Rs < H
    error('̽������С�������ܳ���')
end
% Vp = 7160;
if H > 500e3
    Vp = fun_Vp(H);          % ƽ̨�ٶ�m/s
end
beta = 2*Vp/fr/d_az;    %%���ϵ��;% Brennan����
beta0 = 4*Vp/fr/lambda;
%% �ɵ�����
alpha1 = 0;            %����γ��deg
eta = 90;              %�������deg
%%����������ָ��

Az0 = 90;
El0 = fun_ELAngle(H,R);                  % Elevation (look-down angle) in rad. Flat earth assumption.
%% Thermal Noise Power Computations �ն���������
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 3;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10)-1);        % Effective Receiver Temperature in Kelvin.
Lr = 2.68;                    % System Losses on receive in dB.
Ts = 10^(Lr/10)*Te;           % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts
sigma2 = 1;                   % Normalized Noise Power in Watts.

%% Clutter Patch Geometry computations �Ӳ�·�����νṹ����
% Rs = 130e3;                % (clutter) range of interest in meters. %%б��
% if H >=500e3
%     Rs = Rs;                % (clutter) range of interest in meters. %%б��
% %     Rs = 1300e3;
% end

dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
Re = 6373000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
graze = fun_GrazeAngle(H,R,Rs);             % Grazing angle at the clutter patch in rad (flat earth model).
%%
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.

%% Clutter Covariance Matrix Computations �Ӳ�Э�������


% Calculate Spatial and Doppler Frequencies for k-th clutter patch.
%% ����ģ��
Rsmax = fun_Rsmax(H);                       %���̽��б��m
Rsu=c/(2*fr);                               %�����ģ��б��m
Nk1 = floor((Rsmax-Rs)/Rsu);                %ǰ��ģ����
Nk2 = floor((Rs-H)/Rsu);                    %����ģ����
Rsk1 = Rs + (0:Nk1)*Rsu;                    %ǰ��ģ��б��m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %����ģ��б��m
if isRu == 0
    Rsk = Rs;
else
    Rsk = [Rsk2,Rsk1];%Rs;%                     %ģ��б��m
end

Rk = fun_Rs2R(H,Rsk);                       %ģ���ؾ�m
grazek = fun_GrazeAngle(H,Rk,Rsk);             % Grazing angle at the clutter patch in rad (flat earth model).
elk = fun_ELAngle(H,Rk);                    %��ģ�����븩����deg
% grazek =  fun_GrazeAngle(H,Rk,Rsk)/180*pi;         %��ģ�����������rad
%%%%%%%%%
% Normalized Doppler Frequency of the k-th clutter patch:
if isRotation == 1
    CrabA = fun_CrabAngle(alpha1,eta,H);
    CrabM = fun_CrabMagnitude(alpha1,eta,H);
else
    CrabA = 0;
    CrabM = 1;
end


% Clutter Steering Vector Pre-allocation:
Rc = zeros(Np*N_az*N_el,Np*N_az*N_el); 
for k = 1:Nc %ͬһ���벻ͬ��λ��Ԫ����
%     k
    Rc_t = zeros(Np*N_az*N_el,Np*N_az*N_el);
%     Vc_t = zeros(M*N_az*N_el,1);
%     for i=1:length(elk)  %%%%����ģ��
        c0 = sin(deg2rad(El0))*cos(deg2rad(Az(k)));
        t0 = c0./sin(deg2rad(elk));
        index=find(t0<=1.0001 & t0>=-1.0001);
        az = real(acos(t0(index)));
        omegac = beta0*CrabM.*sin(deg2rad(elk(index))).*cos(az + CrabA);
        for i_az = 1:length(az)
            %% ��������
            At=fun_Ax( N_az, N_el, d_az, d_el, lambda, Az0, El0, rad2deg(az(i_az)), elk(index(i_az)));
            Ar=fun_Ax( N_az, N_el, d_az, d_el, lambda, Az0, El0, rad2deg(az(i_az)), elk(index(i_az)));
            At=abs(At);
            Ar=abs(Ar);
%             % Calculate the Full Array Transmit Power Gain:
%             Gtgain = 10^(Gt/10)*abs(At).^2;
%             % Calculate the Element Receive Power Gain:
%             grgain = 10^(Gel/10)*abs(Ar).^2;
%             % Clutter Patch RCS Calculation: %%�������ӳ������ 
%             PatchArea = Rsk(index(i_az))*dphi*dR*sec(grazek(index(i_az))/180*pi); %�Ӳ���� ,dphi,Ӧ���ǲ�����Ȱ�
%             sigma0 = gamma*sin(grazek(index(i_az))/180*pi);
%             sigma = sigma0*PatchArea;
%             ksi = Pt*Gtgain.*grgain*10^(Gr/10)*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rsk(index(i_az))^4);
           %% ����ʸ�����
            a = fun_RectSteer(N_az, N_el, d_az, d_el, lambda ,rad2deg(az(i_az)), elk(index(i_az)))/sqrt(N_az*N_el);% Spatial Steering Vector.   
            b = exp(1j*pi*omegac(i_az)*(0:Np-1)).'/sqrt(Np); % Temporal Steering Vector
            Vc_t = (At*Ar)*kron(b,a);           % Space-Time Steering Vector.
            Rc_t = Rc_t + (Vc_t * Vc_t'); 
        end
%     end
    Rc = Rc + Rc_t;
end
% Rc = Vc*Vc';                            % Eq. (64)
Rn = eye(Np*N_az*N_el);%sigma2*
% Rn = zeros(Np*N_az*N_el);%sigma2*
Rcn = Rc + Rn;
end

