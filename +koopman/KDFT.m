function [lambdas, Modes, Amps] = KDFT(Snapshots, dt)
%KDFT Compute Koopman modes using Discrete Fourier Transform (FFT), as used
%by Mezic group at UC Santa Barbara.
%
% [lambdas, Modes, Amps] = KDFT( Snapshots, dt )
%    Compute Koopman modes of data in Snapshots matrix. Columns of Snapshots are
%    measurements taken dt apart.
%
%    Koopman modes are taken by computing FFT of every *row* of Snapshots
%    matrix, and then sorting those frequencies that give contribution to a
%    large number of observables.
%
%    lambdas -- list of complex Koopman frequencies, real part is the
%    decay rate, imaginary part (angular) frequency.
%    Modes  -- each L2-normalized column of the matrix is a Koopman mode, corresponding
%    to the lambdas at the same index.
%    Amps   -- optimal L2 amplitudes used to sort the modes in descending order

% Copyright 2015 under BSD license (see LICENSE file).

import koopman.*

% number of snapshots
N = size(Snapshots, 2);

% sampling frequency
Fs = 1/dt;

% Evaluate FFT and rescale it appropriately
F = fft( Snapshots, [], 2 );
Np = size(F,2);
F = F/Np;

% use single-sided FFT
Modes = F(:,1:(floor(Np/2)+1));
Modes(:, 2:end-1) = 2*Modes(:, 2:end-1);

% Frequencies at which modes are obtained
lambdas = complex(0, 2*pi*Fs*(0:(Np/2))/Np).';

%%
% Sort modes according to their optimal L2 contributions
[~,lambdas, Modes, Amps] = sortmodes( lambdas, Modes, Snapshots, dt );
