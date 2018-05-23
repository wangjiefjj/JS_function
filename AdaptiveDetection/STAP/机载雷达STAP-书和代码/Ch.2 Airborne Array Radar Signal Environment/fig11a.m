%% Figure 11a. Illustrating The Clutter Signal Fourier Power Spectrum.

%%
% Coded by Ilias Konsoulas, 3 Aug. 2016.
% Code provided for educational purposes only. All rights reserved.

clc; clear; close all;

%% Radar System Operational Parameters.
fo = 450e6;                   % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4e6;                     % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
M = 18;                       % Number of Pulses per CPI:
Tp = 200e-6;                  % Pulse Width in sec.
N = 18;                       % Number of Array Antenna Elements
Gel = 4;                      % Element Gain in dB
be = -30;                     % Element Backlobe Level in db
Nc = 361;                     % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % Interelement Spacing

% Azimuth angle in degrees:
phi = -180:180;
Lphi = length(phi);
f = zeros(1,Lphi);
AF = zeros(1,Lphi);           % Array Factor vector pre-allocation.

%% Platform Parameters.
beta = 1;                     % Beta Parameter.
ha = 9e3;                     % Platform altitude in meters.
Rc = 13e4;                    % (clutter) range of interest in meters.

%% Thermal Noise Power Computations.
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F = 3;                        % Receiver Noise Figure in dB;
Te = To*(10^(F/10) - 1);      % Effective Receiver Temperature in Kelvin.
Nn = k*Te;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts
sigma2 = 1;                   % Normalized Noise Power.

%% Clutter Patch Geometry Computations.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                   % Radar Range Resolution in meters.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
theta = psi;
phia = 0;                     % Velocity Misalignment angle in degrees.

%% Calculate the Voltage Element Pattern:
for i =1:Lphi
    if abs(phi(i))<=90
        f(i) = cos(phi(i)*pi/180);
    else
        f(i) = 10^(be/10)*cos(phi(i)*pi/180);
    end
end

%% Calculate and Plot the Array Factor (AF) (Voltage):
steering_angle = 0; % Angle of beam steering in degrees.
for k=1:Lphi
    AF(k) = sum(exp(-1i*2*pi/lambda*d*(0:N-1)*(sin(phi(k)*pi/180) ...
                - sin(steering_angle*pi/180))).*cos(phi(k)*pi/180));
end

%% Calculate and Plot the Full Array Transmit Power Gain:
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate and Plot the Element Receive Power Gain:
grgain = 10^(Gel/10)*10^(Gr/10)*abs(f).^2;

%% Clutter Patch RCS Calculation:
PatchArea = Rc*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

%% Calculate and Plot the Clutter to Noise Ration (CNR) for each clutter patch:
ksi = Pt*Gtgain.*grgain*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rc^4);  % Eq. (58)

%% Create Spatial Steering Vector:
a = zeros(N,Nc);
b = zeros(M,Nc);
Vc = zeros(M*N,Nc);

Ksic = sigma2*diag(ksi);

% Platform Velocity for beta parameter value:
va = round(beta*d*fr/2);
Ita = d/lambda*cos(theta);

% Calculate Spatial and Doppler Frequencies.
% Spatial frequency of the k-th clutter patch.
fsp = Ita*sin(phi*pi/180);
% Normalized Doppler Frequency:
omegac = beta*Ita*sin(phi*pi/180 + phia*pi/180);

%% Compute the Clutter Covariance Matrix.

for k=1:Nc
    a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1));                % Spatial Steering Vector.
    b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1));             % Temporal Steering Vector
    Vc(:,k) = kron(b(:,k),a(:,k));                       % Space-Time Steering Vector.
end

Rc = Vc*Ksic*Vc';                                        % Eq. (64)

%% Compute the Clutter Power Spectrum.

phi = -90:90;     Lphi = length(phi);
fd = -150:150;  Lfd = length(fd);
fsp = d/lambda*cos(theta)*sin(phi*pi/180);
omegac = fd/fr;
Pw1 = zeros(Lfd,Lphi);
for m=1:Lphi
    for n=1:Lfd
        a = exp(1i*2*pi*fsp(m)*(0:N-1));                 % Dummy Spatial Steering Vector.
        b = exp(1i*2*pi*omegac(n)*(0:M-1));              % Dummy Doppler Steering Vector
        v = kron(b,a).';
        Pw1(n,m) = (v'*Rc*v)/(M*N);
    end
end

%% Normalization.
max_value = max(max(Pw1));
Pw = Pw1/max_value;

%% Cropping Extreme Values.
[rows cols] = find(10*log10(abs(Pw))<-80);
for i=1:length(rows)
    Pw(rows(i),cols(i)) = 10^(-80/10);
end

%% Plot the Clutter Power Spectrum.
figure('NumberTitle', 'off','Name', 'Figure 11a. Clutter Power Spectrum - Side Looking Radar', ... 
          'Position', [50 50 1150 500]);
subplot(1,2,1);
[Az Doppler] = meshgrid(sin(phi*pi/180),fd);
colormap jet;
pcolor(Az, Doppler/fr, 10*log10(abs(Pw)));
shading interp;
subplot(1,2,2);
[Az Doppler] = meshgrid(sin(phi*pi/180),fd);
colormap jet;
surfc(Az, Doppler/fr, 10*log10(abs(Pw)));
shading interp;
xlabel('sin(Azimuth)');
ylabel('Doppler Frequency (Hz)');
h = colorbar;
% h = colorbar('YTickLabel',{-80:10:0});
set(get(h,'YLabel'),'String','Relative Power (dB)');
tightfig;


