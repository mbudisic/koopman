function [lambda, Modes] = KDFT(Snapshots, dt)
%KDFT Compute Koopman modes using Discrete Fourier Transform (FFT).
%

% number of snapshots
N = size(Snapshots, 2);

% sampling frequency
Fs = 1/dt;

% Evaluate FFT and rescale it appropriately
F = fft( Snapshots, [], 2 );
Np = size(F,2);
F = F/Np;

% use single-sided FFT
Modes = F(:,1:(Np/2+1));
Modes(:, 2:end-1) = 2*Modes(:, 2:end-1);

% Frequencies at which modes are obtained
lambda = complex(0, 2*pi*Fs*(0:(Np/2))/Np).';

% sort modes
[lambda, Modes] = sortmodes( lambda, Modes );
