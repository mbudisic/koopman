function [lambdas, Modes, Amps] = KDFT(Snapshots, dt)
%KDFT Compute Koopman modes using Discrete Fourier Transform (FFT).
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
%    Amps   -- optimal L2 amplitudes used to sort the modes in descending
%    order
%
% The function returns both conjugate pairs of the modes and their
% frequencies.
%
% Note: this method currently correctly operates only on non-transient
% data.
%
% See also DMD_DUKE, DMD, L2OPTIMALMODEAMPLITUDES

% Copyright 2015 under BSD license (see LICENSE file).

import koopman.*

% number of snapshots
N = size(Snapshots, 2);

% sampling frequency
Fs = 1/dt;

% Evaluate FFT and rescale it appropriately
F = fft( Snapshots, N, 2 );
Np = size(F,2);
F = F/Np;


endIdx = (Np-1)/2; % maximum mode index
validateattributes(endIdx,{'numeric'},{'positive','integer'});
assert( size(F,2) == (endIdx*2+1), 'Mode index incorrectly computed');

% use double-sided FFT
idx = [0:endIdx, -endIdx:-1];
lambdas = complex(0, 2*pi*Fs*(idx/Np)).';
Modes = F;

%%
% Compute optimal L2 reconstruction amplitudes
% and sort modes according to their modulus
[~,lambdas, Modes, Amps] = sortmodes( lambdas, Modes, Snapshots, dt );
