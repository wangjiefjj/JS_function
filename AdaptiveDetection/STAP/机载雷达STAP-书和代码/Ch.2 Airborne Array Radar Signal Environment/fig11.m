%% Figure 11. Illustrating Brennan's Rule: Clutter Eigenspectra eigenspectra for the example radar system with different platform velocities.

%%
% Coded by Ilias Konsoulas, 16 Dec. 2014.
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
beta = [0.6 1 2 2.83 3];      % Beta Parameter Vector.
ha = 9e3;                     % Platform altitude in meters.
Rc = 13e4;                    % (clutter) range of interest in meters.

%% Thermal Noise Power Computations.
k = 1.3806488e-23;            % Boltzmann Constant in J/K.
To = 290;                     % Standard room Temperature in Kelvin.
F   = 3;                      % Receiver Noise Figure in dB;
Te = To*(10^(F/10) - 1);      % Effective Receiver Temperature in Kelvin.
Nn = k*Te;                    % Receiver Noise PSD in Watts/Hz.
Pn = Nn*B;                    % Receiver Noise Power in Watts

%% Clutter Patch Geometry Computations.
dphi = 2*pi/Nc;               % Azimuth angle increment in rad.
dR = c/2/B;                   % Radar Range Resolution in meters.
Re = 6370000;                 % Earth Radius in meters.
ae = 4/3*Re;                  % Effective Earth Radius in meters.
psi = asin(ha/Rc);            % Grazing angle at the clutter patch in rad (flat earth model).
gamma = 10^(-3/10);           % Terrain-dependent reflectivity factor.
theta = psi;

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
PatchArea = Rc*dphi*dR*sec(psi);
sigma0 = gamma*sin(psi);
sigma = sigma0*PatchArea;

%% Calculate and Plot the Clutter to Noise Ration (CNR) for each clutter patch:
ksi = Pt*Gtgain.*grgain*lambda^2*sigma/((4*pi)^3*Pn*10^(Ls/10)*Rc^4);

%% Create Spatial Steering Vector:
a = zeros(N,Nc);
b = zeros(M,Nc);
Vc = zeros(M*N,Nc);
Rc = zeros(M*N,M*N);
Ksic = diag(ksi);
colors = [0 0 1; 0 1 0; 1 0 0 ; 1 1 0; 0 1 1;];

figure('NumberTitle', 'off','Name', ...
    'Figure 11. Illustrating Brennan''s Rule: Clutter Eigenspectra for Different Platform Velocities',...
    'Position', [50 50 700 550]);

for i=1:length(beta)
    for k=1:Nc
        fsp = d/lambda*cos(theta)*sin(phi(k)*pi/180); % Spatial frequency of the k-th clutter patch.
        a(:,k) = exp(1i*2*pi*fsp*(0:N-1));            % Spatial Steering Vector.
        omegac = beta(i)*fsp;                         % Normalized Doppler frequency of the k-th clutter patch.
        b(:,k) = exp(1i*2*pi*omegac*(0:M-1));         % Temporal Steering Vector
        Vc(:,k) = kron(b(:,k),a(:,k));                % Space-Time Steering Vector.
    end
    
    Rc = Vc*Ksic*Vc';
    
    plot(10*log10(abs(eig(Rc))),'--s','LineWidth',1,'Color', colors(i,:), ...
        'MarkerEdgeColor','k','MarkerFaceColor',colors(i,:), 'MarkerSize',5);
    hold on;
end

va = round(beta*d*fr/2);

legend(['\beta = 0.6,  v_a = ',num2str(va(1))], ['\beta = 1,     v_a = ',num2str(va(2))], ...
       ['\beta = 2,    v_a = ',num2str(va(3))], ['\beta = 2.83, v_a = ',num2str(va(4))], ...
       ['\beta = 3,    v_a = ',num2str(va(5))]);

ylim([-60 80]); xlim([1 100]);
grid on;
xlabel('Eigenvalue Number');
ylabel('Relative Power (dB)');

for i=1:length(beta)
    X = [round(N+(M-1)*beta(i)), round(N+(M-1)*beta(i))];
    Y = [-60, 80];
    line(X,Y,'Color',colors(i,:),'LineWidth',2)
end
