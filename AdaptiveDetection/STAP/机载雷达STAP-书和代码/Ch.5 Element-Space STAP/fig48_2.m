%% Figure 48.2. Adapted Pattern for Adjacent-Bin post-Doppler STAP for K=2.
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
% Jamming Covariance Matrix:
Rj = kron(eye(M),Phi_j);                                  % Eq. (45)

%% Total Interference Covariance Matrix
Ru = Rc + Rj + Rn;                                        % Eq. (98)
InvRu = inv(Ru);

%% Target Space-Time Steering Vector
phit = 0; thetat = 0;                              % Target azimuth and elevation angles in degrees.
fdt = 100;                                              % Target Doppler Frequency.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);
omegat = fdt/fr;
bt = exp(-1i*2*pi*omegat*(0:M-1)).';  % Target Doppler Steering Vector.
at = exp(-1i*2*pi*fspt*(0:N-1)).';          % Target Spatial Steering Vector.
ta = chebwin(N,30);                              % 30 dB Chebychev Spatial Tapper.
gt = kron(bt,ta.*at);

%% Doppler Filter Bank Creation:
dopplerfilterbank = linspace(-150,150,M+1);
omegadopplerbank = dopplerfilterbank/fr;

%% Doppler Filter Matrix Construction for Adjacent Bin Post-Doppler method:
K = 2;
P = floor(K/2);

U2 = zeros(M,M);
if  mod(K,2)    % If K is odd
    for m=1:M
        U2(:,m) =  1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1));     % Doppler Filter Steering Vector
    end
else            % if K is even:
    outomegadoppler = zeros(1,M);
    for m=1:M
        outomegadoppler(m) = (omegadopplerbank(m) + omegadopplerbank(m+1))/2;
        U2(:,m) =   1/sqrt(M)*exp(-1i*2*pi*outomegadoppler(m)*(0:M-1));     % Doppler Filter Steering Vector
    end
end

td0ab   = ones(M,1);                                                   % Uniform Doppler Taper.
td30ab = chebwin(M,30);                                               % 30-dB Chebyshev Doppler Taper.
td60ab = chebwin(M,60);                                               % 60-dB Chebyshev Doppler Taper.
td90ab = chebwin(M,90);                                               % 90-dB Chebyshev Doppler Taper.

% Create Doppler Filter Bank in Fab matrix for Adjacent Bin post Doppler method:
Fab0   = diag(td0ab)*U2;                           % Eq. 230 without complex conjugation operator.
Fab30 = diag(td30ab)*U2;
Fab60 = diag(td60ab)*U2;
Fab90 = diag(td90ab)*U2;

Fmab0   = zeros(M,K,M);
Fmab30 = zeros(M,K,M);
Fmab60 = zeros(M,K,M);
Fmab90 = zeros(M,K,M);

for m=1:M
    if mod(K,2)                  % if K is odd and >1.
        if (m-P>0) && (m+P<=M)
            %                 [m m-P:m+P]
            %                 omegadopplerbank(m)
            %                omegadopplerbank(m-P:m+P)
            Fmab0(:,:,m)   = Fab0(:,m-P:m+P);                                           % Eq. 231.
            Fmab30(:,:,m) = Fab30(:,m-P:m+P);
            Fmab60(:,:,m) = Fab60(:,m-P:m+P);
            Fmab90(:,:,m) = Fab90(:,m-P:m+P);
        elseif (m-P<=0) && (m+P<=M)
            %                 [m M+(m-P):M  1:m+P]
            %                 omegadopplerbank(m)
            %                 omegadopplerbank([M+(m-P):M  1:m+P])
            Fmab0(:,:,m)   = [Fab0(:,M+(m-P):M)      Fab0(:,1:m+P)];                   % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,M+(m-P):M) Fab30(:,1:m+P)];
            Fmab60(:,:,m) = [Fab60(:,M+(m-P):M) Fab60(:,1:m+P)];
            Fmab90(:,:,m) = [Fab90(:,M+(m-P):M) Fab90(:,1:m+P)];
        elseif m+P>M
            %                  [m m-P:M  1:m+P-M]
            %                  omegadopplerbank(m)
            %                  omegadopplerbank([m-P:M  1:m+P-M])
            Fmab0(:,:,m)   = [Fab0(:,m-P:M)      Fab0(:,1:m+P-M)];                        % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,m-P:M) Fab30(:,1:m+P-M)];
            Fmab60(:,:,m) = [Fab60(:,m-P:M) Fab60(:,1:m+P-M)];
            Fmab90(:,:,m) = [Fab90(:,m-P:M) Fab90(:,1:m+P-M)];
        end
        
    else      % if K is even.
        
        if (m-P>0) && (m+P<=M+1)
            %                     [m m-P:m+P-1]
            %                     omegadopplerbank(m)
            %                     outomegadoppler(m-P:m+P-1)
            Fmab0(:,:,m)   = Fab0(:,m-P:m+P-1);                                          % Eq. 231.
            Fmab30(:,:,m) = Fab30(:,m-P:m+P-1);
            Fmab60(:,:,m) = Fab60(:,m-P:m+P-1);
            Fmab90(:,:,m) = Fab90(:,m-P:m+P-1);
        elseif (m-P<=0) && (m+P<=M)
            %                     [m M+(m-P):M  1:m+P-1]
            %                     omegadopplerbank(m)
            %                     outomegadoppler([M+(m-P):M  1:m+P-1])
            Fmab0(:,:,m)   = [Fab0(:,M+(m-P):M)   Fab0(:,1:m+P-1)];              % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,M+(m-P):M) Fab30(:,1:m+P-1)];
            Fmab60(:,:,m) = [Fab60(:,M+(m-P):M) Fab60(:,1:m+P-1)];
            Fmab90(:,:,m) = [Fab90(:,M+(m-P):M) Fab90(:,1:m+P-1)];
        elseif m+P>M+1
            %                     [m m-P:M 1:m-M+P-1]
            %                     omegadopplerbank(m)
            %                     outomegadoppler([m-P:M 1:m-M+P-1])
            Fmab0(:,:,m)   = [Fab0(:,m-P:M)     Fab0(:,1:m-M+P-1)];           % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,m-P:M)   Fab30(:,1:m-M+P-1)];
            Fmab60(:,:,m) = [Fab60(:,m-P:M)   Fab60(:,1:m-M+P-1)];
            Fmab90(:,:,m) = [Fab90(:,m-P:M)   Fab90(:,1:m-M+P-1)];
        end
    end
    
end
%%     Adjacent-Bin SINR Computations:
m = 16;                % This is the Target's Doppler Bin.

f0abm   = Fmab0(:,:,m);                % Select a Cluster of K Adjacent Doppler Bins.
f30abm = Fmab30(:,:,m);
f60abm = Fmab60(:,:,m);
f90abm = Fmab90(:,:,m);

R0abum   = kron(f0abm,eye(N))'*Ru*kron(f0abm,eye(N));
R30abum = kron(f30abm,eye(N))'*Ru*kron(f30abm,eye(N));
R60abum = kron(f60abm,eye(N))'*Ru*kron(f60abm,eye(N));
R90abum = kron(f90abm,eye(N))'*Ru*kron(f90abm,eye(N));

gt0abm   =  kron(f0abm,eye(N))'*gt;                                             % Eq. 141
gt30abm =  kron(f30abm,eye(N))'*gt;
gt60abm =  kron(f60abm,eye(N))'*gt;
gt90abm =  kron(f90abm,eye(N))'*gt;

% Calculate K*N X 1 Adaptive Weight for m-th Doppler Bin.
w0abm   = R0abum\gt0abm;                                                        % Eq. 205
w30abm = R30abum\gt30abm;
w60abm = R60abum\gt60abm;
w90abm = R90abum\gt90abm;

wab0   = kron(f0abm,eye(N))*w0abm;
wab30 = kron(f30abm,eye(N))*w30abm;
wab60 = kron(f60abm,eye(N))*w60abm;
wab90 = kron(f90abm,eye(N))*w90abm;

%% Adapted Patterns
phi = -90:90; Lphi = length(phi);
fd = -150:150;   Lfd = length(fd);
fsp = d/lambda*cos(theta*pi/180)*sin(phi*pi/180);
omega = fd/fr;
Pw0 = zeros(Lfd,Lphi);
Pw30 = zeros(Lfd,Lphi);
Pw60 = zeros(Lfd,Lphi);
Pw90 = zeros(Lfd,Lphi);
for m1=1:Lphi
    for n=1:Lfd
        a = exp(-1i*2*pi*fsp(m1)*(0:N-1));           % Dummy Spatial Steering Vector.
        b = exp(-1i*2*pi*omega(n)*(0:M-1));     % Dummy Doppler Steering Vector
        v = kron(b,a).';
        Pw0(n,m1)   = abs(wab0'*v)^2;
        Pw30(n,m1) = abs(wab30'*v)^2;
        Pw60(n,m1) = abs(wab60'*v)^2;
        Pw90(n,m1) = abs(wab90'*v)^2;
    end
end

%% Normalisation:
max_value0   = max(max(Pw0));
max_value30 = max(max(Pw30));
max_value60 = max(max(Pw60));
max_value90 = max(max(Pw90));

Pw0 = Pw0/max_value0;
Pw30 = Pw30/max_value30;
Pw60 = Pw60/max_value60;
Pw90 = Pw90/max_value90;

[rows0 cols0] = find(10*log10(abs(Pw0))<-150);
for i=1:length(rows0)
    Pw0(rows0(i),cols0(i)) = 10^(-150/10);
end

[rows30 cols30] = find(10*log10(abs(Pw30))<-150);
for i=1:length(rows30)
    Pw30(rows30(i),cols30(i)) = 10^(-150/10);
end

[rows60 cols60] = find(10*log10(abs(Pw60))<-150);
for i=1:length(rows60)
    Pw60(rows60(i),cols60(i)) = 10^(-150/10);
end

[rows90 cols90] = find(10*log10(abs(Pw90))<-150);
for i=1:length(rows90)
    Pw90(rows90(i),cols90(i)) = 10^(-150/10);
end
%% Plot the Adapted Pattern:
figure('NumberTitle', 'off','Name', ...
    ['Figure 48.2. Adapted Patterns for Adjacent-Bin post Doppler STAP for K = ', ...
    num2str(K), ' and Doppler bin ' num2str(m)], ...
    'Position', [1 1 1000 1000]);
subplot(2,2,1);
[Az Doppler] = meshgrid(sin(phi*pi/180),fd);
colormap jet;
pcolor(Az, Doppler, 10*log10(abs(Pw0)));
shading interp;
xlim([-1 1])
ylim([-150 150]);
% xlabel('sin(Azimuth)');
ylabel('Doppler Frequency (Hz)');
h = colorbar;
% set(get(h,'YLabel'),'String','Relative Power (dB)');
title('Doppler Filters Untapered');

subplot(2,2,2);
pcolor(Az, Doppler, 10*log10(abs(Pw30)));
shading interp;
xlim([-1 1])
ylim([-150 150]);
% xlabel('sin(Azimuth)');
% ylabel('Doppler Frequency (Hz)');
title('30 dB Chebychev Doppler Taper');
h = colorbar;
set(get(h,'YLabel'),'String','Relative Power (dB)');

subplot(2,2,3);
pcolor(Az, Doppler, 10*log10(abs(Pw60)));
shading interp;
xlim([-1 1])
ylim([-150 150]);
xlabel('sin(Azimuth)');
ylabel('Doppler Frequency (Hz)');
title('60 dB Chebychev Doppler Taper');
h = colorbar;
% set(get(h,'YLabel'),'String','Relative Power (dB)');

subplot(2,2,4);
pcolor(Az, Doppler, 10*log10(abs(Pw90)));
shading interp;
xlim([-1 1])
ylim([-150 150]);
xlabel('sin(Azimuth)');
% ylabel('Doppler Frequency (Hz)');
title('90 dB Chebychev Doppler Taper');
h = colorbar;
set(get(h,'YLabel'),'String','Relative Power (dB)');

tightfig;
