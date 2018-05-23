%% Figure 15. Clutter Eigenspectra for Different Misalignment Angles.
%%
%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
% Code provided for educational purposes only. All rights reserved.

clc; clear; close all;

%% Radar System Operational Parameters.
fo = 450e6;                   % Operating Frequency in Hz
Pt = 200e3;                   % Peak Transmit Power 200 kW
Gt = 22;                      % Transmit Gain in dB
Gr = 10;                      % Column Receive Gain in dB
B  = 4*1e6;                   % Receiver Instantaneous Bandwidth in Hz
Ls = 4;                       % System Losses in dB
fr = 300;                     % PRF in Hz
M = 18;                       % Number of Pulses per CPI.
Tp = 200*1e-6;                % Pulse Width in sec.
N = 18;                       % Number Array Antenna Elements
Gel = 4;                      % Element Gain in dB.

% Transmit Taper: Uniform
be = -30;                     % Element Backlobe Level in db
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
Pn = Nn*B;                    % Receiver Noise Power in Watts

%% Clutter Patch Geometry Parameters.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                   % Radar Range Resolution in meters.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
psi = asin(ha/Rcik);          % Grazing angle at the clutter patch in rad (flat earth model).
theta = psi;                  % Elevation (look-down angle). Flat earth assumption.
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
phia = [0 1 10 45 90];        % Velocity Misalignment angle in degrees.

colors = [0 0 1; 0 1 0; 1 0 0 ; 1 1 0; 0 1 1; 1 0 1;];

%% Calculate the Voltage Element Pattern.
for i =1:Lphi
    if abs(phi(i))<=90
        f(i) = cos(phi(i)*pi/180);
    else
        f(i) = 10^(be/10)*cos(phi(i)*pi/180);
    end
end

%% Calculate and Plot the Array Factor (AF) (Voltage).
steering_angle = 0; % Angle of beam steering in degrees.
for k=1:Lphi
    AF(k) = sum(exp(-1i*2*pi/lambda*d*(0:N-1)*(sin(phi(k)*pi/180) ...
        - sin(steering_angle*pi/180))).*cos(phi(k)*pi/180));
end

%% Calculate and Plot the Full Array Transmit Power Gain.
Gtgain = 10^(Gt/10)*abs(AF).^2;

% Calculate and Plot the Element Receive Power Gain:
grgain = 10^(Gel/10)*10^(Gr/10)*abs(f).^2;

%% Clutter Patch RCS Calculation.
PatchArea = Rcik*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

%% Calculate and Plot the Clutter to Noise Ration (CNR) for each clutter patch.
ksi = Pt*Gtgain.*grgain*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rcik^4);
Ksic = diag(ksi);

% Platform Velocity for beta parameter value:
va = beta*d*fr/2;
Ita = d/lambda*cos(theta);

figure('NumberTitle', 'off','Name', ...
    'Figure 15. Clutter Eigenspectra for Different Misalignment Angles','Position',[1 1 700 600]);

for i=1:length(phia)
    % Steering Vector Pre-allocation:
    a = zeros(N,Nc);
    b = zeros(M,Nc);
    vc = zeros(M*N,Nc);
    Vc = zeros(M*N,M*N);

    % Spatial frequency of the k-th clutter patch.
    fsp = Ita*sin(phi*pi/180);                                 % Eq. (83a)

    % Normalized Doppler Frequency:
    omegac = beta*Ita*sin(phi*pi/180 + phia(i)*pi/180);        % Eq. (83b)
    
    for k=1:Nc
        % Normalized Doppler Frequency fold-over       .
        if omegac(k)>0.5 && omegac(k)<=1.5
            omegac(k) = omegac(k) - 1;
        elseif omegac(k)<-0.5 && omegac(k)>=-1.5
            omegac(k) = omegac(k) + 1;
        end
        
        if omegac(k)>1.5 && omegac(k)<=2.5
            omegac(k) = omegac(k) - 2;
        elseif omegac(k)<-1.5 && omegac(k)>=-2.5
            omegac(k) = omegac(k) + 2;
        end
        
        if omegac(k)>2.5 && omegac(k)<=3.5
            omegac(k) = omegac(k) - 3;
        elseif omegac(k)<-2.5 && omegac(k)>=-3.5
            omegac(k) = omegac(k) + 3;
        end
        
        a(:,k) = exp(1i*2*pi*fsp(k)*(0:N-1)); % Spatial Steering Vector.
        b(:,k) = exp(1i*2*pi*omegac(k)*(0:M-1)); % Temporal Steering Vector
        Vc(:,k) = kron(b(:,k),a(:,k));  % Space-TIme Steering Vector.
    end
    
    Rc = Vc*Ksic*Vc';                                          % Eq. (64)
    
    plot(10*log10(abs(eig(Rc))),'--s','LineWidth',1,'Color', colors(i,:), ...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(i,:), 'MarkerSize',5);
    hold on;
end
ylim([-80 80]); xlim([0 100]);
grid on;
legend('\phi_a = 0 deg', '\phi_a = 1 deg', '\phi_a = 10 deg', '\phi_a = 45 deg', '\phi_a = 90 deg');
ylabel('Relative Power (dB)'); xlabel('Eigenvalue Number');

