%% Figure 48. SINR Loss performance for PRI-staggered and adjacent-bin post-Doppler STAP.
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

%% Target Space-Time Steering Vector:
phit = 0; thetat = 0;                                     % Target azimuth and elevation angles in degrees.
fdt = 100;                                                % Target Doppler Frequency.
fspt = d/lambda*cos(thetat*pi/180)*sin(phit*pi/180);
omegat = fdt/fr;
at = exp(-1i*2*pi*fspt*(0:N-1)).';                        % Target Spatial Steering Vector.
ta = chebwin(N,30);                                       % 30 dB Chebychev Spatial Tapper.

%% Doppler Filter Bank Creation:
dopplerfilterbank = linspace(0,300,M+1);
omegadopplerbank = dopplerfilterbank/fr;
fd = 0:.5:300;   Lfd = length(fd);
omegad = fd/fr;
SNRo = M*N;

%% Doppler Filter Matrix Construction for PRI-Staggered Post-Doppler method
K = 2;
P = floor(K/2);
M1= M - K +1;

U1 = zeros(M1,M);
for m=1:M
    U1(:,m) =  1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M1-1));  % Doppler Filter Steering Vector
end

td0   = ones(M1,1);
td30 = chebwin(M1,30);                                                % 30-dB Chebyshev Doppler Taper.
td60 = chebwin(M1,60);                                                % 60-dB Chebyshev Doppler Taper.
td90 = chebwin(M1,90);                                                % 90-dB Chebyshev Doppler Taper.

F0   = diag(td0)*U1;                                                  % Eq. 227.
F30 = diag(td30)*U1;
F60 = diag(td60)*U1;
F90 = diag(td90)*U1;

%% Create Doppler Filter Bank in Fm Matrix for PRI-Staggered Post-Doppler method:
Fm0  = zeros(M,K,M);
Fm30 = zeros(M,K,M);
Fm60 = zeros(M,K,M);
Fm90 = zeros(M,K,M);
for m=1:M
    Fm0(:,:,m)  = toeplitz([F0(:,m);  zeros(K-1,1)],[F0(1,m)  zeros(1,K-1)]);              % Eq. 229.
    Fm30(:,:,m) = toeplitz([F30(:,m); zeros(K-1,1)],[F30(1,m) zeros(1,K-1)]);
    Fm60(:,:,m) = toeplitz([F60(:,m); zeros(K-1,1)],[F60(1,m) zeros(1,K-1)]);
    Fm90(:,:,m) = toeplitz([F90(:,m); zeros(K-1,1)],[F90(1,m) zeros(1,K-1)]);
end

%% Doppler Filter Matrix Construction for Adjacent Bin Post-Doppler method
U2 = zeros(M,M);
if  mod(K,2)    % If K is odd
    for m=1:M
        U2(:,m) =  1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1));  % Doppler Filter Steering Vector
    end
else            % if K is even:
    outomegadoppler = zeros(1,M);
    for m=1:M
        outomegadoppler(m) = (omegadopplerbank(m) + omegadopplerbank(m+1))/2;
        U2(:,m) =   1/sqrt(M)*exp(-1i*2*pi*outomegadoppler(m)*(0:M-1));  % Doppler Filter Steering Vector
    end
end

td0ab   = ones(M,1);                                                  % Uniform Doppler Taper.
td30ab = chebwin(M,30);                                               % 30-dB Chebyshev Doppler Taper.
td60ab = chebwin(M,60);                                               % 60-dB Chebyshev Doppler Taper.
td90ab = chebwin(M,90);                                               % 90-dB Chebyshev Doppler Taper.

% Create Doppler Filter Bank in Fab matrix for Adjacent Bin post Doppler method:
Fab0   = diag(td0ab)*U2;                                % Eq. 230 without complex conjugation operator.
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
            Fmab0(:,:,m)   = Fab0(:,m-P:m+P);                                               % Eq. 231.
            Fmab30(:,:,m) = Fab30(:,m-P:m+P);
            Fmab60(:,:,m) = Fab60(:,m-P:m+P);
            Fmab90(:,:,m) = Fab90(:,m-P:m+P);
        elseif (m-P<=0) && (m+P<=M)
            %                 [m M+(m-P):M  1:m+P]
            %                 omegadopplerbank(m)
            %                 omegadopplerbank([M+(m-P):M  1:m+P])
            Fmab0(:,:,m)   = [Fab0(:,M+(m-P):M)      Fab0(:,1:m+P)];                      % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,M+(m-P):M) Fab30(:,1:m+P)];
            Fmab60(:,:,m) = [Fab60(:,M+(m-P):M) Fab60(:,1:m+P)];
            Fmab90(:,:,m) = [Fab90(:,M+(m-P):M) Fab90(:,1:m+P)];
        elseif m+P>M
            %                  [m m-P:M  1:m+P-M]
            %                  omegadopplerbank(m)
            %                  omegadopplerbank([m-P:M  1:m+P-M])
            Fmab0(:,:,m)   = [Fab0(:,m-P:M)      Fab0(:,1:m+P-M)];                          % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,m-P:M) Fab30(:,1:m+P-M)];
            Fmab60(:,:,m) = [Fab60(:,m-P:M) Fab60(:,1:m+P-M)];
            Fmab90(:,:,m) = [Fab90(:,m-P:M) Fab90(:,1:m+P-M)];
        end
        
    else      % if K is even.
        
        if (m-P>0) && (m+P<=M+1)
            %                     [m m-P:m+P-1]
            %                     omegadopplerbank(m)
            %                     outomegadoppler(m-P:m+P-1)
            Fmab0(:,:,m)   = Fab0(:,m-P:m+P-1);                                      % Eq. 231.
            Fmab30(:,:,m) = Fab30(:,m-P:m+P-1);
            Fmab60(:,:,m) = Fab60(:,m-P:m+P-1);
            Fmab90(:,:,m) = Fab90(:,m-P:m+P-1);
        elseif (m-P<=0) && (m+P<=M)
            %                     [m M+(m-P):M  1:m+P-1]
            %                     omegadopplerbank(m)
            %                     outomegadoppler([M+(m-P):M  1:m+P-1])
            Fmab0(:,:,m)   = [Fab0(:,M+(m-P):M)      Fab0(:,1:m+P-1)];                % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,M+(m-P):M) Fab30(:,1:m+P-1)];
            Fmab60(:,:,m) = [Fab60(:,M+(m-P):M) Fab60(:,1:m+P-1)];
            Fmab90(:,:,m) = [Fab90(:,M+(m-P):M) Fab90(:,1:m+P-1)];
        elseif m+P>M+1
            %                     [m m-P:M 1:m-M+P-1]
            %                     omegadopplerbank(m)
            %                     outomegadoppler([m-P:M 1:m-M+P-1])
            Fmab0(:,:,m)   = [Fab0(:,m-P:M)        Fab0(:,1:m-M+P-1)];           % Eq. 231.
            Fmab30(:,:,m) = [Fab30(:,m-P:M)   Fab30(:,1:m-M+P-1)];
            Fmab60(:,:,m) = [Fab60(:,m-P:M)   Fab60(:,1:m-M+P-1)];
            Fmab90(:,:,m) = [Fab90(:,m-P:M)   Fab90(:,1:m-M+P-1)];
        end
    end
    
end

%% LSINR Computation for Optimum Fully Adaptive Case
LSINRopt = zeros(1,Lfd);
for n=1:Lfd
    bt = exp(-1i*2*pi*omegad(n)*(0:M-1)).'; % Target Doppler Steering Vector
    vt = kron(bt,at);
    w = InvRu*vt; %#ok<MINV>
    LSINRopt(n) = real(w'*vt)/SNRo;
end

%% LSINR Computations for 0, 30, 60, and 90 dB Chebyshev Taper.
SINR0_mat    = zeros(M,Lfd);
SINR30_mat   = zeros(M,Lfd);
SINR60_mat   = zeros(M,Lfd);
SINR90_mat   = zeros(M,Lfd);
SINRab0_mat  = zeros(M,Lfd);
SINRab30_mat = zeros(M,Lfd);
SINRab60_mat = zeros(M,Lfd);
SINRab90_mat = zeros(M,Lfd);

for m=1:M
    %% PRI-Staggered SINR Computations
    f0m  = Fm0(:,:,m);
    f30m = Fm30(:,:,m);
    f60m = Fm60(:,:,m);
    f90m = Fm90(:,:,m);
    
    R0um  = kron(f0m,eye(N))'*Ru*kron(f0m,eye(N));
    R30um = kron(f30m,eye(N))'*Ru*kron(f30m,eye(N));
    R60um = kron(f60m,eye(N))'*Ru*kron(f60m,eye(N));
    R90um = kron(f90m,eye(N))'*Ru*kron(f90m,eye(N));
    
    bdfb = exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1)).';
    gt = kron(bdfb,ta.*at);                   % Desired Vector common for both methods.
    
    gt0m  = kron(f0m,eye(N))'*gt;
    gt30m = kron(f30m,eye(N))'*gt;
    gt60m = kron(f60m,eye(N))'*gt;
    gt90m = kron(f90m,eye(N))'*gt;
    
    w0m  = R0um\gt0m;                      % Calculate K*N X 1 Adaptive Weight for m-th Doppler Bin.
    w30m = R30um\gt30m;
    w60m = R60um\gt60m;
    w90m = R90um\gt90m;
    
    w0  = kron(f0m,eye(N))*w0m;
    w30 = kron(f30m,eye(N))*w30m;
    w60 = kron(f60m,eye(N))*w60m;
    w90 = kron(f90m,eye(N))*w90m;
    
    for n=1:Lfd
        bt = exp(1i*2*pi*omegad(n)*(0:M-1)).'; % Dummy Target Doppler Steering Vector
        vt = kron(bt,at);
        SINR0_mat(m,n)  = abs(w0'*vt)^2/real(w0'*Ru*w0);
        SINR30_mat(m,n) = abs(w30'*vt)^2/real(w30'*Ru*w30);
        SINR60_mat(m,n) = abs(w60'*vt)^2/real(w60'*Ru*w60);
        SINR90_mat(m,n) = abs(w90'*vt)^2/real(w90'*Ru*w90);
    end
    
    %% Adjacent-Bin SINR Computations
    f0abm  = Fmab0(:,:,m);                              % Select a Cluster of K Adjacent Doppler Bins.
    f30abm = Fmab30(:,:,m);
    f60abm = Fmab60(:,:,m);
    f90abm = Fmab90(:,:,m);
    
    R0abum  = kron(f0abm,eye(N))'*Ru*kron(f0abm,eye(N));
    R30abum = kron(f30abm,eye(N))'*Ru*kron(f30abm,eye(N));
    R60abum = kron(f60abm,eye(N))'*Ru*kron(f60abm,eye(N));
    R90abum = kron(f90abm,eye(N))'*Ru*kron(f90abm,eye(N));
    
    gt0abm  = kron(f0abm,eye(N))'*gt;                                             % Eq. 141
    gt30abm = kron(f30abm,eye(N))'*gt;
    gt60abm = kron(f60abm,eye(N))'*gt;
    gt90abm = kron(f90abm,eye(N))'*gt;
    
    % Calculate K*N X 1 Adaptive Weight for m-th Doppler Bin.
    w0abm  = R0abum\gt0abm;                                                        % Eq. 205
    w30abm = R30abum\gt30abm;
    w60abm = R60abum\gt60abm;
    w90abm = R90abum\gt90abm;
    
    wab0   = kron(f0abm,eye(N))*w0abm;
    wab30 = kron(f30abm,eye(N))*w30abm;
    wab60 = kron(f60abm,eye(N))*w60abm;
    wab90 = kron(f90abm,eye(N))*w90abm;
    
    for n=1:Lfd
        bt = exp(1i*2*pi*omegad(n)*(0:M-1)).'; % Dummy Target Doppler Steering Vector
        vt = kron(bt,at);
        SINRab0_mat(m,n)  = abs(wab0'*vt)^2/(wab0'*Ru*wab0);
        SINRab30_mat(m,n) = abs(wab30'*vt)^2/(wab30'*Ru*wab30);
        SINRab60_mat(m,n) = abs(wab60'*vt)^2/(wab60'*Ru*wab60);
        SINRab90_mat(m,n) = abs(wab90'*vt)^2/(wab90'*Ru*wab90);
    end
    
end

%% PRI-Staggered SINR Loss Computations
SINR0   = max(abs(SINR0_mat));
SINR30 = max(abs(SINR30_mat));
SINR60 = max(abs(SINR60_mat));
SINR90 = max(abs(SINR90_mat));
LSINR0   = SINR0/SNRo;
LSINR30 = SINR30/SNRo;
LSINR60 = SINR60/SNRo;
LSINR90 = SINR90/SNRo;

%% Adjacent Bin SINR Loss Computations
SINRab0  = max(abs(SINRab0_mat));
SINRab30 = max(abs(SINRab30_mat));
SINRab60 = max(abs(SINRab60_mat));
SINRab90 = max(abs(SINRab90_mat));
LSINRab0  = SINRab0/SNRo;
LSINRab30 = SINRab30/SNRo;
LSINRab60 = SINRab60/SNRo;
LSINRab90 = SINRab90/SNRo;

%% Plot the SINR Losses
figure('NumberTitle', 'off','Name', ...
       'Figure 48. SINR loss performance for PRI-Staggered and Adjacent-Bin post-Doppler STAP',...
       'Position',[1 1 1000 1000]);
subplot(2,2,1);
plot(fd,10*log10(LSINRopt),'LineWidth',1.5)
hold on;
plot(fd,10*log10(LSINR0),'r','LineWidth',1.5)
plot(fd,10*log10(LSINRab0),'g','LineWidth',1.5)
title('Untapered (uniform) Doppler Filters');
ylabel('SINR Loss (dB)');
xlabel('Target Doppler Frequency (Hz)');
ylim([-30 1]);
xlim([-5 305]);
hleg1 = legend('Optimum', ['PRI-Staggered K=',num2str(K)], ['Adjacent Bin K=',num2str(K)], ...
               'Location', 'South');
set(hleg1,'FontSize',8);
grid on;

subplot(2,2,2);
plot(fd,10*log10(LSINRopt),'LineWidth',1.5)
hold on;
plot(fd,10*log10(LSINR30),'r','LineWidth',1.5)
plot(fd,10*log10(LSINRab30),'g','LineWidth',1.5)
title('30-dB Chebyshev Doppler Filters');
% ylabel('SINR Loss (dB)');
xlabel('Target Doppler Frequency (Hz)');
ylim([-30 1]);
xlim([-5 305]);
hleg2 = legend('Optimum', ['PRI-Staggered K=',num2str(K)], ['Adjacent Bin K=',num2str(K)],...
               'Location', 'South');
set(hleg2,'FontSize',8);
grid on;

subplot(2,2,3);
plot(fd,10*log10(LSINRopt),'LineWidth',1.5)
hold on;
plot(fd,10*log10(LSINR60),'r','LineWidth',1.5)
plot(fd,10*log10(LSINRab60),'g','LineWidth',1.5)
title('60-dB Chebyshev Doppler Filters');
ylabel('SINR Loss (dB)');
xlabel('Target Doppler Frequency (Hz)');
ylim([-30 1]);
xlim([-5 305]);
hleg3 = legend('Optimum', ['PRI-Staggered K=',num2str(K)], ['Adjacent Bin K=',num2str(K)],...
                'Location', 'South');
set(hleg3,'FontSize',8);
grid on;

subplot(2,2,4);
plot(fd,10*log10(LSINRopt),'LineWidth',1.5)
hold on;
plot(fd,10*log10(LSINR90),'r','LineWidth',1.5)
plot(fd,10*log10(LSINRab90),'g','LineWidth',1.5)
title('90-dB Chebyshev Doppler Filters');
% ylabel('SINR Loss (dB)');
xlabel('Target Doppler Frequency (Hz)');
ylim([-30 1]);
xlim([-5 305]);
hleg4 = legend('Optimum', ['PRI-Staggered K=',num2str(K)], ['Adjacent Bin K=',num2str(K)],...
                'Location', 'South');
set(hleg4,'FontSize',8);
grid on;

tightfig;