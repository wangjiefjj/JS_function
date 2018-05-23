%% Figure 47. Clutter Eigenspectra for multiwindow post-Doppler approaches with K = 2.

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

%% Doppler Filter Matrix Construction for PRI-Staggered Post-Doppler method:
dopplerfilterbank = linspace(0,300,M+1);
omegadopplerbank = dopplerfilterbank/fr;
K = 2;
P = floor(K/2);
M1= M - K +1;

U1 = zeros(M1,M);
for m=1:M
    U1(:,m) = 1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M1-1));
end

td0   = ones(M1,1);
td30 = chebwin(M1,30);                                               % 30-dB Chebyshev Doppler Taper.
td60 = chebwin(M1,60);                                               % 60-dB Chebyshev Doppler Taper.
td90 = chebwin(M1,90);                                               % 90-dB Chebyshev Doppler Taper.
F0    = diag(td0)*U1;
F30 = diag(td30)*U1;
F60 = diag(td60)*U1;
F90 = diag(td90)*U1;

%% Solve M Separate N-dimensional Adaptive Problems for PRI-Staggered Post-Doppler:
Rcm0   = zeros(N*K,N*K,M);
Rcm30 = zeros(N*K,N*K,M);
Rcm60 = zeros(N*K,N*K,M);
Rcm90 = zeros(N*K,N*K,M);
for m=1:M
    Fm0   = toeplitz([F0(:,m);   zeros(K-1,1)],[F0(1,m) zeros(1,K-1)]);                  % Eq. 229.
    Fm30 = toeplitz([F30(:,m); zeros(K-1,1)],[F30(1,m) zeros(1,K-1)]);
    Fm60 = toeplitz([F60(:,m); zeros(K-1,1)],[F60(1,m) zeros(1,K-1)]);
    Fm90 = toeplitz([F90(:,m); zeros(K-1,1)],[F90(1,m) zeros(1,K-1)]);
    
    Rcm0(:,:,m)   = kron(Fm0  ,eye(N))'*Rc*kron(Fm0,   eye(N));                         % Eq. 214.
    Rcm30(:,:,m) = kron(Fm30,eye(N))'*Rc*kron(Fm30,eye(N));
    Rcm60(:,:,m) = kron(Fm60,eye(N))'*Rc*kron(Fm60,eye(N));
    Rcm90(:,:,m) = kron(Fm90,eye(N))'*Rc*kron(Fm90,eye(N));
end


%% Doppler Filter Matrix Construction for Adjacent Bin Post-Doppler method:
U2 = zeros(M,M);
if  mod(K,2)    % If K is odd
    for m=1:M
        U2(:,m) = 1/sqrt(M)*exp(-1i*2*pi*omegadopplerbank(m)*(0:M-1)); % Doppler Filter Steering Vector
    end
else            % while if K is even:
    for m=1:M
        U2(:,m) = 1/sqrt(M)*exp(-1i*2*pi*(omegadopplerbank(m) + omegadopplerbank(m+1))/2*(0:M-1)); 
    end
end

td0   = ones(M,1);
td30 = chebwin(M,30);                                               % 30-dB Chebyshev Doppler Taper.
td60 = chebwin(M,60);                                               % 60-dB Chebyshev Doppler Taper.
td90 = chebwin(M,90);                                               % 90-dB Chebyshev Doppler Taper.
Fab0    = diag(td0)*U2;
Fab30 = diag(td30)*U2;
Fab60 = diag(td60)*U2;
Fab90 = diag(td90)*U2;

%% Solve M Separate N-dimensional Adaptive Problems for Adjacent Bin Post-Doppler:
Rcmab0   = zeros(N*K,N*K,M);
Rcmab30 = zeros(N*K,N*K,M);
Rcmab60 = zeros(N*K,N*K,M);
Rcmab90 = zeros(N*K,N*K,M);

for m=1:M
    if mod(K,2)                  % if K is odd.
        if (m-P>0) && (m+P<=M)
            Fmab0   = Fab0(:,m-P:m+P);                                        % Eq. 231.
            Fmab30 = Fab30(:,m-P:m+P);
            Fmab60 = Fab60(:,m-P:m+P);
            Fmab90 = Fab90(:,m-P:m+P);
        elseif (m-P<=0) && (m+P<=M)
            Fmab0   = [Fab0(:,M+(m-P):M)  Fab0(:,1:m+P)];                    % Eq. 231.
            Fmab30 = [Fab30(:,M+(m-P):M) Fab30(:,1:m+P)];
            Fmab60 = [Fab60(:,M+(m-P):M) Fab60(:,1:m+P)];
            Fmab90 = [Fab90(:,M+(m-P):M) Fab90(:,1:m+P)];
        elseif m+P>M
            Fmab0   = [Fab0(:,m-P:M)  Fab0(:,1:m+P-M)];                        % Eq. 231.
            Fmab30 = [Fab30(:,m-P:M) Fab30(:,1:m+P-M)];
            Fmab60 = [Fab60(:,m-P:M) Fab60(:,1:m+P-M)];
            Fmab90 = [Fab90(:,m-P:M) Fab90(:,1:m+P-M)];
        end
        
    else      % if K is even.
        
        if (m-P>0) && (m+P<=M+1)
            Fmab0   = Fab0(:,m-P:m+P-1);                                   % Eq. 231.
            Fmab30 = Fab30(:,m-P:m+P-1);
            Fmab60 = Fab60(:,m-P:m+P-1);
            Fmab90 = Fab90(:,m-P:m+P-1);
        elseif (m-P<=0) && (m+P<=M)
            Fmab0   = [Fab0(:,M+(m-P):M)  Fab0(:,1:m+P-1)];                % Eq. 231.
            Fmab30 = [Fab30(:,M+(m-P):M) Fab30(:,1:m+P-1)];
            Fmab60 = [Fab60(:,M+(m-P):M) Fab60(:,1:m+P-1)];
            Fmab90 = [Fab90(:,M+(m-P):M) Fab90(:,1:m+P-1)];
        elseif m+P>M+1
            Fmab0   = [Fab0(:,m-P:M) Fab0(:,1:m-M+P-1)];                    % Eq. 231.
            Fmab30 = [Fab30(:,m-P:M) Fab30(:,1:m-M+P-1)];
            Fmab60 = [Fab60(:,m-P:M) Fab60(:,1:m-M+P-1)];
            Fmab90 = [Fab90(:,m-P:M) Fab90(:,1:m-M+P-1)];
        end
    end
    
    
    Rcmab0(:,:,m)   = kron(Fmab0  ,eye(N))'*Rc*kron(Fmab0, eye(N));        % Eq. 214.
    Rcmab30(:,:,m) = kron(Fmab30,eye(N))'*Rc*kron(Fmab30,eye(N));
    Rcmab60(:,:,m) = kron(Fmab60,eye(N))'*Rc*kron(Fmab60,eye(N));
    Rcmab90(:,:,m) = kron(Fmab90,eye(N))'*Rc*kron(Fmab90,eye(N));
    
end

%% Plot the Clutter Eigenspectra for PRI-Staggered Post-Doppler for various bins:
figure('NumberTitle', 'off','Name', ...
    ' Figure 47. Clutter Eigenspectra for multiwindow post-Doppler Approaches with K=2',...
    'Position',[1 1 900 1200]);
subplot(3,2,1);
bins = [1 4 10];
plot(10*log10(sort(abs(eig(Rcm0(:,:,bins(1)))),'descend')),'b.-')
hold on;
plot(10*log10(sort(abs(eig(Rcm30(:,:,bins(1)))),'descend')),'r.-')
plot(10*log10(sort(abs(eig(Rcm60(:,:,bins(1)))),'descend')),'g.-')
plot(10*log10(sort(abs(eig(Rcm90(:,:,bins(1)))),'descend')),'c.-')
title('PRI-Staggered, Whitened, Bin #0');
ylabel('Relative Power (dB)');
ylim([-80 80]);
xlim([1 36]);
hleg1 = legend('Uniform','30 dB','60 dB','90 dB');
set(hleg1,'FontSize',8);
grid on;

subplot(3,2,3);
plot(10*log10(sort(abs(eig(Rcm0(:,:,bins(2)))),'descend')),'b.-')
hold on;
plot(10*log10(sort(abs(eig(Rcm30(:,:,bins(2)))),'descend')),'r.-')
plot(10*log10(sort(abs(eig(Rcm60(:,:,bins(2)))),'descend')),'g.-')
plot(10*log10(sort(abs(eig(Rcm90(:,:,bins(2)))),'descend')),'c.-')
title('PRI-Staggered, Whitened, Bin #3');
ylabel('Relative Power (dB)');
ylim([-80 80]);
xlim([1 36]);
hleg2 = legend('Uniform','30 dB','60 dB','90 dB');
set(hleg2,'FontSize',8);
grid on;

subplot(3,2,5);
plot(10*log10(sort(abs(eig(Rcm0(:,:,bins(3)))),'descend')),'b.-')
hold on;
plot(10*log10(sort(abs(eig(Rcm30(:,:,bins(3)))),'descend')),'r.-')
plot(10*log10(sort(abs(eig(Rcm60(:,:,bins(3)))),'descend')),'g.-')
plot(10*log10(sort(abs(eig(Rcm90(:,:,bins(3)))),'descend')),'c.-')
title('PRI-Staggered, Whitened, Bin #9');
ylabel('Relative Power (dB)');
ylim([-80 80]);
xlim([1 36]);
hleg3 = legend('Uniform','30 dB','60 dB','90 dB');
set(hleg3,'FontSize',8);
grid on;

% Plot the clutter eigenspectra for Adjacent Bin Post-Doppler.
subplot(3,2,2);
plot(10*log10(sort(abs(eig(Rcmab0(:,:,bins(1)))),'descend')),'b.-')
hold on;
plot(10*log10(sort(abs(eig(Rcmab30(:,:,bins(1)))),'descend')),'r.-')
plot(10*log10(sort(abs(eig(Rcmab60(:,:,bins(1)))),'descend')),'g.-')
plot(10*log10(sort(abs(eig(Rcmab90(:,:,bins(1)))),'descend')),'c.-')
title('Adjacent Bin, Whitened, Bin #0');
% ylabel('Relative Power (dB)');
ylim([-80 80]);
xlim([1 36]);
hleg1 = legend('Uniform','30 dB','60 dB','90 dB');
set(hleg1,'FontSize',8);
grid on;

subplot(3,2,4);
plot(10*log10(sort(abs(eig(Rcmab0(:,:,bins(2)))),'descend')),'b.-')
hold on;
plot(10*log10(sort(abs(eig(Rcmab30(:,:,bins(2)))),'descend')),'r.-')
plot(10*log10(sort(abs(eig(Rcmab60(:,:,bins(2)))),'descend')),'g.-')
plot(10*log10(sort(abs(eig(Rcmab90(:,:,bins(2)))),'descend')),'c.-')
title('Adjacent Bin, Whitened, Bin #3');
% ylabel('Relative Power (dB)');
ylim([-80 80]);
xlim([1 36]);
hleg2 = legend('Uniform','30 dB','60 dB','90 dB');
set(hleg2,'FontSize',8);
grid on;

subplot(3,2,6);
plot(10*log10(sort(abs(eig(Rcmab0(:,:,bins(3)))),'descend')),'b.-')
hold on;
plot(10*log10(sort(abs(eig(Rcmab30(:,:,bins(3)))),'descend')),'r.-')
plot(10*log10(sort(abs(eig(Rcmab60(:,:,bins(3)))),'descend')),'g.-')
plot(10*log10(sort(abs(eig(Rcmab90(:,:,bins(3)))),'descend')),'c.-')
title('Adjacent Bin, Whitened, Bin #9');
%ylabel('Relative Power (dB)');
ylim([-80 80]);
xlim([1 36]);
hleg3 = legend('Uniform','30 dB','60 dB','90 dB');
set(hleg3,'FontSize',8);
grid on;

tightfig;