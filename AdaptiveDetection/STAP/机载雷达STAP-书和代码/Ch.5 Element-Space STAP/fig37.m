%% Figure 37. Composite pattern for Doppler bin 5. 
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
K = 2;                        % Number of Pulses per sub-CPI.
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
b = zeros(K,Nc);
Vc = zeros(K*N,Nc);

for k=1:Nc
    a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1));       % Spatial Steering Vector.
    b(:,k) = exp(1i*2*pi*omegac(k)*(0:K-1));    % Temporal Steering Vector
    Vc(:,k) = kron(b(:,k),a(:,k));              % Space-Time Steering Vector.
end

Rcsub = Vc*Ksic*Vc';
Rnsub = sigma2*eye(K*N);

%% Jammer Covariance Matrix Calculation.
J = 2;                                                   % Number of Jammers.
thetaj = 0; phij = [-40 25];                             % Jammer elevation and azimuth angles in degrees.
R_j = [370 370]*1e3;
Sj = 1e-3;                                               % Jammer ERPD in Watts/Hz.
fspj = d/lambda*cos(thetaj*pi/180)*sin(phij*pi/180);     % Spatial frequency of the j-th jammer.
Lrj = 1.92;                                              % System Losses on Receive in dB.
Aj = zeros(N,J);
for j=1:J
    Aj(:,j) =  exp(1i*2*pi*fspj(j)*(0:N-1));             % Jammer Spatial Steering Vector.
end

indices= zeros(1,J);
for j=1:J
    indices(j) = find(phi == phij(j));
end
grgn = grgain(indices);
ksi_j = (Sj*grgn*lambda^2)./((4*pi)^2.*Nn*10^(Lrj/10).*R_j.^2);    

Ksi_j = sigma2*diag(ksi_j);
Phi_j = Aj*Ksi_j*Aj';                                    % Eq. (47)
Rjsub = kron(eye(K),Phi_j);

%% Analytic Interference Covariance Matrix for first sub-CPI:
Rusub = Rcsub + Rjsub + Rnsub;

%% Target Space-Time Steering Vector:
phit = 0; thetat = 0;                                   % Target azimuth and elevation angles in degrees.
fdt = 100;                                              % Target Doppler Frequency.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);    % Target Spatial Frequency
at = exp(1i*2*pi*fspt*(0:N-1)).';                       % Target Spatial Steering Vector.
ta = chebwin(N,30);                                     % 30 dB Chebychev Spatial Tapper.
bt = [1; -1];                                           % K-pulse Temporal Target Steering Vector
tb = [1;  1];                                           % (K x 1) Binomial Taper for K = 2.

gt = kron(tb.*bt,ta.*at);                               % sub-CPI desired response.

%% Tapered, Element-Space STAP Solution
wsub = Rusub\gt;

%% W matrix construction. Equation (180).
M1 = M - K + 1;                                         % Number of sub-CPI's.
W = zeros(M*N,M1);

% Assume that weight vectors for each sub-CPI are equal.
for k=1:M1
     W((k-1)*N+1:(k+1)*N,k) = wsub;
end

td =  chebwin(M1,40);                                   % 40 dB Chebychev Doppler Taper.
fd5 = 93.75;                                            % Center Frequency of Doppler Bin no. 5.
omegad5 = fd5/fr;                                       % Normalized Center Frequency of Doppler Bin no. 5.

um = 1/sqrt(M1)*exp(1i*2*pi*omegad5*(0:M1-1)).';
fm = td.*um;
wm = W*fm;                                              % Eq. (184). m-th Doppler bin composite weight vector.

%% Adapted Patterns
phi = -90:.5:90; Lphi = length(phi);
fd = -150:.5:150;   Lfd = length(fd);
fsp = d/lambda*cos(theta*pi/180)*sin(phi*pi/180);
omega = fd/fr;
Pw1 = zeros(Lfd,Lphi);
Pw2 = zeros(Lfd,Lphi);
for m=1:Lphi
    for n=1:Lfd
        a = exp(1i*2*pi*fsp(m)*(0:N-1));               % Dummy Spatial Steering Vector.
        b = exp(1i*2*pi*omega(n)*(0:M-1));             % Dummy Doppler Steering Vector
        v = kron(b,a).';
        Pw2(n,m) = abs(wm'*v)^2;
    end
end

%% Normalisation
max_value = max(max(Pw2));
Pw = Pw2/max_value;
[rows cols] = find(10*log10(abs(Pw))<-100);
for i=1:length(rows)
    Pw(rows(i),cols(i)) = 10^(-100/10);
end

%% Plot the Adapted Pattern for Doppler Bin #5
figure('NumberTitle', 'off','Name', 'Figure 37. Composite pattern for Doppler bin 5. (a) Full pattern.', ...
       'Position', [1 1 700 600]);
[Az Doppler] = meshgrid(sin(phi*pi/180),fd);
colormap jet;
pcolor(Az, Doppler, 10*log10(abs(Pw)));
shading interp;
xlim([-1 1])
ylim([-150 150]);
xlabel('sin(Azimuth)');
ylabel('Doppler Frequency (Hz)');
h = colorbar;
set(get(h,'YLabel'),'String','Relative Power (dB)');

%% Plot the Principal Cuts:
figure('NumberTitle', 'off','Name', ...
       'Figure 37. Composite Pattern for Doppler Bin 5. (b) Principal Cuts at Target Azimuth and Doppler', ...
       'Position', [1 1 700 600]);
% a. Cut of the Adapted Pattern at Doppler = 100 Hz.
subplot(2,1,1);
plot(sin(phi*pi/180), 10*log10(abs(Pw(fd == 100,:))));
ylim([-100 0.5]); xlim([-1  1]);
ylabel('Relative Power (dB)');
xlabel('sin(Azimuth)');
title('Doppler Frequency = 100 Hz');
grid on;

% b. Cut of the Adapted Pattern at Target Azimuth = 0 degrees.
subplot(2,1,2);
plot(fd, 10*log10(abs(Pw(:,phi == phit))));
ylim([-100 0.5]); xlim([-150 150]);
ylabel('Relative Power (dB)');
xlabel('Target Doppler Frequency (Hz)');
title('Azimuth = 0 deg');
grid on;
