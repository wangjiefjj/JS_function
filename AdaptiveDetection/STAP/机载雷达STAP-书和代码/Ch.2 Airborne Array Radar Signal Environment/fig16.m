%% Figure 16. Clutter Eigenspectra with Misalignment for Different Backlobe Levels.
 
%%
%  Coded by Ilias Konsoulas, 16 Dec. 2014.
%  Code provided for educational purposes only. All rights reserved.

clc; clear; close all;

%% Radar System Operational Parameters:
fo = 450e6;                   % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4*1e6;                   % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
M = 18;                       % Number of Pulses per CPI.
Tp = 200*1e-6;                % Pulse Width in sec.
N = 18;                       % Number of Array Antenna Elements
Gel = 4;                      % Element Gain in dB.
be = -20:-10:-70;             % Element Backlobe Level in db.
Nc = 361;                     % Number of clutter patches uniformly distributed in azimuth.
c   = 299792458;              % Speed of Light in m/sec.
lambda = c/fo;                % Operating wavelength in meters.
d = lambda/2;                 % Interelement Spacing

% Azimuth angle in degrees:
phi = -180:180;
Lphi = length(phi);
f = zeros(1,Lphi);
AF = zeros(1,Lphi);           % Array Factor pre-allocation.

%% Platform Parameters:
beta = 1;                     % beta parameter.
ha = 9e3;                     % Platform altitude in meters.
Rcik = 13e4;                  % (clutter) range of interest in meters.

%% Thermal Noise Power Computations.
k = 1.3806488*1e-23;          % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F   = 3;                      % Receiver Noise Figure in dB;
Te = To*(10^(F/10) - 1);      % Effective Receiver Temperature in Kelvin.
Nn = k*Te;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts.
sigma2 = 1;                   % Normalized Noise Power.

%% Clutter Patch Geometry Parameters.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                   % Radar Range Resolution in meters.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
psi = asin(ha/Rcik);          % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Elevation (look-down angle). Flat earth assumption.
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
phia = 10;                    % Velocity Misalignment angle in degrees.

colors = [0 0 1; 0 1 0; 1 0 0 ; 1 1 0; 0 1 1; 1 0 1;];

%% Calculate the Array Factor (AF) (Voltage).
steering_angle = 0; % Angle of beam steering in degrees.
for k=1:Lphi
    AF(k) = sum(exp(-1i*2*pi/lambda*d*(0:N-1)*(sin(phi(k)*pi/180) ...
                 - sin(steering_angle*pi/180))).*cos(phi(k)*pi/180));
end

%% Calculate the Full Array Transmit Power Gain.
Gtgain = 10^(Gt/10)*abs(AF).^2;

%% Clutter Patch RCS Calculation.
PatchArea = Rcik*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

%% Calculate and Plot the Clutter Eigenspectrum for Different Backlobe Levels.
figure('NumberTitle', 'off','Name', ...
       'Figure 16. Clutter Eigenspectra with Misalignment for Different Backlobe Levels', 'Position',[1 1 750 600]);

for i1=1:length(be)
    
    % Calculate the Voltage Element Pattern:
    for i =1:Lphi
        if abs(phi(i))<=90
            f(i) = cos(phi(i)*pi/180);
        else
            f(i) = 10^(be(i1)/10)*cos(phi(i)*pi/180);
        end
    end
    
    % Calculate the Element Receive Power Gain:
    grgain = 10^(Gel/10)*10^(Gr/10)*abs(f).^2;
    
    ksi = Pt*Gtgain.*grgain*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rcik^4);   % Eq. (58)
    Ksic = sigma2*diag(ksi);                                                  % Eq. (66)
    
    % Platform Velocity for beta parameter value:
    va = beta*d*fr/2;
    Ita = d/lambda*cos(theta);
    
    % Steering Vector Pre-allocation:
    a = zeros(N,Nc);
    b = zeros(M,Nc);
    Vc = zeros(M*N,M*N);
    % Spatial frequency of the k-th clutter patch.
    fsp = Ita*sin(phi*pi/180);
    % Normalized Doppler Frequency:
    omegac = beta*Ita*sin(phi*pi/180 + phia*pi/180);
    
    for k=1:Nc
        a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1)); % Spatial Steering Vector.
        b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1)); % Time Steering Vector
        Vc(:,k) = kron(b(:,k),a(:,k));  % Space-TIme Steering Vector.
    end
    
    Rc = Vc*Ksic*Vc';                                                         % Eq. (64);
    
    plot(10*log10(abs(eig(Rc))),'--s','LineWidth',1,'Color', colors(i1,:), ...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(i1,:), 'MarkerSize',5);
    hold on;
end

ylim([-80 80]); xlim([0 100]);
grid on;
legend('b_e = -20 dB', 'b_e = -30 dB', 'b_e = -40 dB', 'b_e = -50 dB', 'b_e = -60 dB', 'b_e = -70 dB');
ylabel('Relative Power (dB)'); xlabel('Eigenvalue Number');
title('Misalignment Angle: \phi_a = 10\circ');

