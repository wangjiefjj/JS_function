function [ Rcn,El0,d_az,lambda,fr,Np,N_az ] = fun_SBRR( isRu, isRotation,R)
%%%根据P156也表6.2的参数设置
if nargin <3
    R = 1300e3;
end
%% 雷达系统参数
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
N_az = 16;  %16               % Number of Array Antenna Elements, 发射时是384,32*12
N_el = 1;
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d_az = lambda/2*1;              % 间隔
% d_az = lambda/2*1.08;              % 间隔
d_el = lambda/2*1.39;              % 间隔
% Azimuth angle in degrees 方位角:
% Az = linspace(-0,179,Nc);
Az = 0:1:180;
EL = 0:1:90;
Nc = length(Az);
LAz = length(Az);
f = zeros(1,LAz);
AF = ones(1,LAz);           % Array Factor pre-allocation. 天线阵因子
Sys_DOF = Np*N_az*N_el;
dR = c/2/B;                   % Radar Range Resolution in meters.
%% Platform Parameters:
H = 506e3;%9e3 %700e3                    % Platform altitude in meters.
Rs = fun_R2Rs(H,R);
if Rs < H
    error('探测距离过小，不可能出现')
end
% Vp = 7160;
if H > 500e3
    Vp = fun_Vp(H);          % 平台速度m/s
end
beta = 2*Vp/fr/d_az;    %%混叠系数;% Brennan因子
beta0 = 4*Vp/fr/lambda;
%% 可调参数
alpha1 = 0;            %卫星纬度deg
eta = 90;              %卫星倾角deg
%%主波束天线指向

Az0 = 90;
El0 = fun_ELAngle(H,R);                  % Elevation (look-down angle) in rad. Flat earth assumption.
%% Thermal Noise Power Computations 终端噪声计算
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 3;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10)-1);        % Effective Receiver Temperature in Kelvin.
Lr = 2.68;                    % System Losses on receive in dB.
Ts = 10^(Lr/10)*Te;           % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts
sigma2 = 1;                   % Normalized Noise Power in Watts.

%% Clutter Patch Geometry computations 杂波路径几何结构计算
% Rs = 130e3;                % (clutter) range of interest in meters. %%斜距
% if H >=500e3
%     Rs = Rs;                % (clutter) range of interest in meters. %%斜距
% %     Rs = 1300e3;
% end

dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
Re = 6373000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
graze = fun_GrazeAngle(H,R,Rs);             % Grazing angle at the clutter patch in rad (flat earth model).
%%
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.

%% Clutter Covariance Matrix Computations 杂波协方差计算


% Calculate Spatial and Doppler Frequencies for k-th clutter patch.
%% 距离模糊
Rsmax = fun_Rsmax(H);                       %最大探测斜距m
Rsu=c/(2*fr);                               %最大无模糊斜距m
Nk1 = floor((Rsmax-Rs)/Rsu);                %前向模糊数
Nk2 = floor((Rs-H)/Rsu);                    %后向模糊数
Rsk1 = Rs + (0:Nk1)*Rsu;                    %前向模糊斜距m
Rsk2 = Rs - (Nk2:-1:1)*Rsu;                 %后向模糊斜距m
if isRu == 0
    Rsk = Rs;
else
    Rsk = [Rsk2,Rsk1];%Rs;%                     %模糊斜距m
end

Rk = fun_Rs2R(H,Rsk);                       %模糊地距m
grazek = fun_GrazeAngle(H,Rk,Rsk);             % Grazing angle at the clutter patch in rad (flat earth model).
elk = fun_ELAngle(H,Rk);                    %各模糊距离俯仰角deg
% grazek =  fun_GrazeAngle(H,Rk,Rsk)/180*pi;         %各模糊距离掠射角rad
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
for k = 1:Nc %同一距离不同方位单元叠加
%     k
    Rc_t = zeros(Np*N_az*N_el,Np*N_az*N_el);
%     Vc_t = zeros(M*N_az*N_el,1);
%     for i=1:length(elk)  %%%%距离模糊
        c0 = sin(deg2rad(El0))*cos(deg2rad(Az(k)));
        t0 = c0./sin(deg2rad(elk));
        index=find(t0<=1.0001 & t0>=-1.0001);
        az = real(acos(t0(index)));
        omegac = beta0*CrabM.*sin(deg2rad(elk(index))).*cos(az + CrabA);
        for i_az = 1:length(az)
            %% 天线增益
            At=fun_Ax( N_az, N_el, d_az, d_el, lambda, Az0, El0, rad2deg(az(i_az)), elk(index(i_az)));
            Ar=fun_Ax( N_az, N_el, d_az, d_el, lambda, Az0, El0, rad2deg(az(i_az)), elk(index(i_az)));
            At=abs(At);
            Ar=abs(Ar);
%             % Calculate the Full Array Transmit Power Gain:
%             Gtgain = 10^(Gt/10)*abs(At).^2;
%             % Calculate the Element Receive Power Gain:
%             grgain = 10^(Gel/10)*abs(Ar).^2;
%             % Clutter Patch RCS Calculation: %%反射因子乘以面积 
%             PatchArea = Rsk(index(i_az))*dphi*dR*sec(grazek(index(i_az))/180*pi); %杂波面积 ,dphi,应该是波束宽度吧
%             sigma0 = gamma*sin(grazek(index(i_az))/180*pi);
%             sigma = sigma0*PatchArea;
%             ksi = Pt*Gtgain.*grgain*10^(Gr/10)*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rsk(index(i_az))^4);
           %% 导向矢量求解
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

