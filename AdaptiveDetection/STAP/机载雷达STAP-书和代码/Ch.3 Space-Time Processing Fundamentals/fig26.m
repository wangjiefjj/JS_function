%% Figure 26. SINR for the optimum and tapered fully adaptive STAP, including Doppler straddling losses.

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

% Clutter Steering Vector Pre-allocation:
a = zeros(N,Nc);
b = zeros(M,Nc);
Vc = zeros(M*N,Nc);

for k=1:Nc
    a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1));    % Spatial Steering Vector.
    b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1)); % Temporal Steering Vector
    Vc(:,k) = kron(b(:,k),a(:,k));           % Space-Time Steering Vector.
end

Rc = Vc*Ksic*Vc';                            % Eq. (64)

Rn = sigma2*eye(M*N);

%% Jamming Covariance Matrix Calculation
J = 2;                                                       % Number of Jammers.
thetaj = 0; phij = [-40 25];                                 % Jammer elevation and azimuth angles in degrees.
R_j = [370 370]*1e3;
Sj = 1e-3;                                                   % Jammer ERPD in Watts/Hz.
fspj = d/lambda*cos(thetaj*pi/180)*sin(phij*pi/180);         % Spatial frequency of the j-th jammer.
Lrj = 1.92;                                                  % System Losses on Receive in dB.
Aj = zeros(N,J);
for j=1:J
    Aj(:,j) =  exp(1i*2*pi*fspj(j)*(0:N-1));                 % Jammer Spatial Steering Vector.
end

indices= zeros(1,J);
for j=1:J
    indices(j) = find(phi == phij(j));
end
grgn = grgain(indices);
ksi_j = (Sj*grgn*lambda^2)./((4*pi)^2.*Nn*10^(Lrj/10).*R_j.^2);

Ksi_j = sigma2*diag(ksi_j);
Phi_j = Aj*Ksi_j*Aj';                                        % Eq. (47)
% Jamming Covariance Matrix:
Rj = kron(eye(M),Phi_j);                                     % Eq. (45)

%% Total Interference Covariance Matrix
Ru = Rc + Rj + Rn;                                           % Eq. (98)

%% SINR calculation for Optimum and Tapered Fully Adaptive STAP:
ta = chebwin(N,30);                                          % 30 dB Chebychev Spatial Tapper.
tb = chebwin(M,40);                                          % 40 dB Chebychev Doppler Taper.
phit = 0; thetat = 0;                                        % Target Azimuth and Elevation Angles.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);         % Target Spatial Frequency.
a = exp(1i*2*pi*fspt*(0:N-1)).';                             % Target Spatial Steering Vector.
fd = 0:.5:300;   Lfd = length(fd);
omega = fd/fr;
dopplerfilterbank = linspace(0,300,M+1);                     % Doppler Filterbank frequencies.
omegadopplerbank = dopplerfilterbank/fr;
SINRopt_mat = zeros(length(dopplerfilterbank),Lfd);
SINRtap_mat = zeros(length(dopplerfilterbank),Lfd);
InvRu = inv(Ru);

for m=1:length(dopplerfilterbank)
     bm = exp(1i*2*pi*omegadopplerbank(m)*(0:M-1)).';        % Doppler Filter Steering Vector
     vm = kron(bm,a);
     gt = kron(tb.*bm,ta.*a);
     wm = InvRu*vm; %#ok<*MINV>                              % Eq. (116)
     wmtap = InvRu*gt;
    
     for n=1:Lfd
         b = exp(1i*2*pi*omega(n)*(0:M-1)).';                % Dummy Target Doppler Steering Vector
         v = kron(b,a);
         SINRopt_mat(m,n) = wm'*v;                           % Eq. (114) 
         SINRtap_mat(m,n) = abs(wmtap'*v)^2/real(wmtap'*gt); % Eq. (115)
     end
end

SINRopt = max(abs(SINRopt_mat));                             % Eq. (117) for Optimum Fully Adaptive Case
SINRtap = max(abs(SINRtap_mat));                             % Eq. (117) for Tapered Fully Adaptive Case

%% Plot the Adapted Pattern:
figure('NumberTitle', 'off','Name', ...
       'Figure 26. SINR for the optimum and tapered fully adaptive STAPs, including Doppler straddling losses.',...
       'Position', [1 1 600 500]);
plot(fd,10*log10(SINRopt),'LineWidth',1.5)
hold on;
plot(fd,10*log10(SINRtap),'r','LineWidth',1.5)
ylabel('SINR (dB)');
xlabel('Target Doppler Frequency (Hz)');
ylim([-5 26]);
xlim([-5 305]);
legend('Optimum Fully Adaptive', 'Tapered Fully Adaptive', 'Location','South')
grid on;