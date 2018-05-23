%% Figure 38. SINR loss for Element Space pre-Doppler. Zero intrinsic cluller motion.
%%
%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
% Code provided for educational purposes only. All rights reserved.

clc; clear; close all;

%% Radar System Operational Parameters
fo = 450e6;                   % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4e6;                     % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
Tr = 1/fr;                    % PRI in sec.
M = 18;                       % Number of Pulses per CPI.
K1 = 2;                       % Number of Pulses per sub-CPI.
K2 = 3;                       % Number of Pulses per sub-CPI.
M1 = M - K1 + 1;              % Number of sub-CPI's for K=2.
M2 = M - K2 + 1;              % Number of sub-CPI's for K=3.
Tp = 200e-6;                  % Pulse Width in sec.
N = 18;                       % Number of Array Antenna Elements
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
Nc = 360;                     % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % Interelement Spacing

% Azimuth angle in degrees:
phi = -180:179;
Lphi = length(phi);
f = zeros(1,Lphi);
AF = zeros(1,Lphi);           % Array Factor pre-allocation.

% Platform Parameters:
beta = 1;                     % beta parameter.
ha = 9e3;                     % Platform altitude in meters.

%% Thermal Noise Power Computations
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 3;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10)-1);        % Effective Receiver Temperature in Kelvin.
Lr = 2.68;                    % System Losses on receive in dB.
Ts = 10^(Lr/10)*Te;           % Reception System Noise Temperature in Kelvin.
Nn = k*Ts;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts
sigma2 = 1;                   % Normalized Noise Power in Watts.

%% Clutter Patch Geometry computations
Rcik = 130000;                % (clutter) range of interest in meters.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                   % Radar Range Resolution in meters.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
psi = asin(ha/Rcik);          % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Elevation (look-down angle) in rad. Flat earth assumption.
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
phia = 0;                     % Velocity Misalignment angle in degrees.

%% Clutter-to-Noise Ratio (CNR) Calculation

% Calculate the Voltage Element Pattern:
for i =1:Lphi
    if abs(phi(i))<=90
        f(i) = cos(phi(i)*pi/180);
    else
        f(i) = 10^(be/10)*cos(phi(i)*pi/180);
    end
end

% Calculate the Array Factor (AF) (Voltage):
steering_angle = 0; % Angle of beam steering in degrees.
for k=1:Lphi
    AF(k) = sum(exp(-1i*2*pi/lambda*d*(0:N-1)*(sin(phi(k)*pi/180) ...
                          - sin(steering_angle*pi/180))));
end

% Calculate the Full Array Transmit Power Gain:
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate the Element Receive Power Gain:
grgain = 10^(Gel/10)*abs(f).^2;

% Clutter Patch RCS Calculation:
PatchArea = Rcik*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

% Calculate the Clutter to Noise Ratio (CNR) for each clutter patch:
ksi = Pt*Gtgain.*grgain*10^(Gr/10)*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rcik^4);
Ksic = sigma2*diag(ksi);

%% Clutter Covariance Matrix Computations

% Platform Velocity for beta parameter value:
va = round(beta*d*fr/2);
Ita = d/lambda*cos(theta);

% Calculate Spatial and Doppler Frequencies for k-th clutter patch.
% Spatial frequency of the k-th clutter patch:
fsp = Ita*sin(phi*pi/180);
% Normalized Doppler Frequency of the k-th clutter patch:
omegac = beta*Ita*sin(phi*pi/180 + phia*pi/180);

% Clutter Steering Vector Pre-allocation for sub-CPI of size K:
a = zeros(N,Nc);
b = zeros(M,Nc);
b1 = zeros(K1,Nc);
b2 = zeros(K2,Nc);
Vc   = zeros(M*N,Nc);
Vc1 = zeros(K1*N,Nc);
Vc2 = zeros(K2*N,Nc);

for k=1:Nc
    a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1));                 % Spatial Steering Vector.
    b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1));              % Temporal Steering Vector.
    b1(:,k) = exp(1i*2*pi*omegac(k)*(0:K1-1));            % Temporal Steering Vector for K = 2.
    b2(:,k) = exp(1i*2*pi*omegac(k)*(0:K2-1));            % Temporal Steering Vector for K = 3.
    Vc(:,k) = kron(b(:,k),a(:,k));                        % Space-Time Steering Vector.
    Vc1(:,k) = kron(b1(:,k),a(:,k));                      % Space-Time Steering Vector for K = 2.
    Vc2(:,k) = kron(b2(:,k),a(:,k));                      % Space-Time Steering Vector for K = 3.
end

Rc = Vc*Ksic*Vc';
Rcsub1 = Vc1*Ksic*Vc1';
Rcsub2 = Vc2*Ksic*Vc2';

Rn = sigma2*eye(M*N);
Rnsub1 = sigma2*eye(K1*N);
Rnsub2 = sigma2*eye(K2*N);

%% Jammer Covariance Matrix Calculation
J = 2;                                                    % Number of Jammers.
thetaj = 0; phij = [-40 25];                              % Jammer elevation and azimuth angles in degrees.
R_j = [370 370]*1e3;
Sj = 1e-3;                                                % Jammer ERPD in Watts/Hz.
fspj = d/lambda*cos(thetaj*pi/180)*sin(phij*pi/180);      % Spatial frequency of the j-th jammer.
Lrj = 1.92;                                               % System Losses on Receive in dB.
Aj = zeros(N,J);
for j=1:J
    Aj(:,j) =  exp(1i*2*pi*fspj(j)*(0:N-1));              % Jammer Spatial Steering Vector.
end

indices= zeros(1,J);
for j=1:J
    indices(j) = find(phi == phij(j));
end
grgn = grgain(indices);
ksi_j = (Sj*grgn*lambda^2)./((4*pi)^2.*Nn*10^(Lrj/10).*R_j.^2);    

Ksi_j = sigma2*diag(ksi_j);
Phi_j = Aj*Ksi_j*Aj';                                     % Eq. (47)

Rj = kron(eye(M),Phi_j);
Rjsub1 = kron(eye(K1),Phi_j);
Rjsub2 = kron(eye(K2),Phi_j);

%% Analytic Interference Covariance Matrix for Fully Adaptive STAP:
Ru = Rc + Rj + Rn;

%% Analytic Interference Covariance Matrix for first sub-CPI (K=2):
Rusub1 = Rcsub1 + Rjsub1 + Rnsub1;

%% Analytic Interference Covariance Matrix for second sub-CPI (K=3):
Rusub2 = Rcsub2 + Rjsub2 + Rnsub2;

%% Target Space-Time Steering Vector:
phit = 0; thetat = 0;                                     % Target azimuth and elevation angles in degrees.
fdt = 100;                                                % Target Doppler Frequency.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);
omegact = fdt/fr;
at = exp(1i*2*pi*fspt*(0:N-1)).';                         % Target Spatial Steering Vector.
ta = chebwin(N,30);                                       % 30 dB Chebychev Spatial Tapper.

bt1 = [1; -1;];
bt2 = [1; -1; 1;];
tb1 = [1;  1;];                                           % Binomial Doppler Taper for K = 2.
tb2 = [1; 2; 1;];                                         % Binomial Doppler Taper for K = 3.

gt1 = kron(tb1.*bt1,ta.*at);                              % sub-CPI #1 desired response.
gt2 = kron(tb2.*bt2,ta.*at);                              % sub-CPI #2 desired response.
 
%% Tapered, Element-Space STAP Solution
wsub1 = Rusub1\gt1;                                       % Weight Vector for K = 2;
wsub2 = Rusub2\gt2;                                       % Weight Vector for K = 3;

%% W1 matrix construction. Equation (180).
W1 = zeros(M*N,M1);
% Assume that weight vectors for each sub-CPI are equal.
for k=1:M1
    W1((k-1)*N+1:(k+1)*N,k) = wsub1;
end

W2 = zeros(M*N,M2);
% W2 matrix construction. Equation (180).
% Assume that weight vectors for each sub-CPI are equal,.
for k=1:M2
    W2((k-1)*N+1:(k+2)*N,k) = wsub2;
end

td1 =  chebwin(M1,40);                                    % 40 dB Chebychev Doppler Taper.
td2 =  chebwin(M2,40);                                    % 40 dB Chebychev Doppler Taper.

%% SINR calculation for Optimum and pre-Doppler Element Space STAP:
phit = 0; thetat = 0;                                     % Target azimuth and elevation angles in degrees.
fdt = 100;                                                % Target Doppler Frequency.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);      % Target Spatial Frequency
at = exp(1i*2*pi*fspt*(0:N-1)).';                         % Target Spatial Steering Vector.
fd = 0:.5:300;   Lfd = length(fd);
omegad = fd/fr;
dopplerfilterbank1 = linspace(0,300,M1+1);
dopplerfilterbank2 = linspace(0,300,M2+1);
omegadopplerbank1 = dopplerfilterbank1/fr;
omegadopplerbank2 = dopplerfilterbank2/fr;
LSINRopt = zeros(1,Lfd);
SINRsub1_mat = zeros(length(dopplerfilterbank1),Lfd);
SINRsub2_mat = zeros(length(dopplerfilterbank2),Lfd);

InvRu = inv(Ru);
SNRo = M*N;

%% LSINR Computation for Optimum Fully Adaptive Case:
for n=1:Lfd
    bt = exp(1i*2*pi*omegad(n)*(0:M-1)).';                      % Dummy Target Doppler Steering Vector
    vt = kron(bt,at);
    w = InvRu*vt; %#ok<MINV>
    LSINRopt(n) = real(w'*vt)/SNRo;
end

%% SINR Computations for sub-CPI with K=2.
for m=1:length(dopplerfilterbank1)
    um = exp(1i*2*pi*omegadopplerbank1(m)*(0:M1-1)).';         % Doppler Filter Steering Vector
    fm1 = td1.*um;
    wm1 = W1*fm1;                                              % m-th Doppler bin composite weight vector.
    
    for n=1:Lfd
        bt = exp(1i*2*pi*omegad(n)*(0:M-1)).';                 % Dummy Target Doppler Steering Vector
        vt = kron(bt,at);
        SINRsub1_mat(m,n) = abs(wm1'*vt)^2/real(wm1'*Ru*wm1);  % Eq. (185)
    end
end

SINRsub1 = max(abs(SINRsub1_mat));                             % Eq. (186)
LSINRsub1 = SINRsub1/SNRo;

%% SINR Computations for sub-CPI with K=3.
for m=1:length(dopplerfilterbank2)
    um = exp(1i*2*pi*omegadopplerbank2(m)*(0:M2-1)).';         % Doppler Filter Steering Vector
    fm2 = td2.*um;
    wm2 = W2*fm2;                                              % m-th Doppler bin composite weight vector.
    
    for n=1:Lfd
        bt = exp(1i*2*pi*omegad(n)*(0:M-1)).';                 % Dummy Target Doppler Steering Vector
        vt = kron(bt,at);
        SINRsub2_mat(m,n) = abs(wm2'*vt)^2/real(wm2'*Ru*wm2);  % Eq. (185)
    end
end

SINRsub2 = max(abs(SINRsub2_mat));                            % Eq. (186)
LSINRsub2 = SINRsub2/SNRo;

%% Plot the SINR Losses
figure('NumberTitle', 'off','Name', ...
       'Figure 38. SINR loss for element space pre-Doppler. Zero intrinsic clutter motion',...
       'Position', [1 1 700 600]);
plot(fd,10*log10(LSINRopt),'LineWidth',1.5)
hold on;
plot(fd,10*log10(LSINRsub1),'r','LineWidth',1.5)
hold on;
plot(fd,10*log10(LSINRsub2),'g','LineWidth',1.5)
ylabel('SINR Loss (dB)');
xlabel('Target Doppler Frequency (Hz)');
ylim([-30 1]);
xlim([-5 305]);
legend('Optimum Fully Adaptive', 'Pre-Doppler K=2','Pre-Doppler K=3', 'Location','South')
grid on;