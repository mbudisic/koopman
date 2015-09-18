function [lambdas, Modes, Amps] = DMD_Duke(Snapshots, dt, varargin)
%DMD_DUKE Compute Koopman modes by Dynamic Mode Decomposition algorithm by Duke et al.
%
% This is the algorithm by Duke, Daniel, Julio Soria, and Damon
% Honnery. 2012. “An Error Analysis of the Dynamic Mode Decomposition.”
% Experiments in Fluids 52 (2): 529–42. doi:10.1007/s00348-011-1235-7.
%
% The algorithm uses the QR decomposition of the input snapshot matrix to
% construct the DMD matrix.
%
% [lambdas, Modes, Amps] = DMD_Duke( Snapshots, dt )
%    Compute DMD of data in Snapshots matrix. Columns of Snapshots are
%    measurements taken dt apart.
%
%    lambdas -- list of complex Dynamic Mode frequencies, real part is the
%    decay rate, imaginary part (angular) frequency.
%    Modes  -- each L2-normalized column of the matrix is a Dynamic Mode, corresponding
%    to the lambdas at the same index.
%    Amps   -- optimal L2 amplitudes used to sort the modes in descending order
%
% [lambdas, Modes, Amps] = DMD_Duke( ..., db ) If db set to true, first
%    "de-bias" using Hemati, Rowley procedure. Alternatively, if db is a
%    positive integer, use db modes to debias.
%
% The function returns both conjugate pairs of the modes and their frequencies.
%
% See also DMD, KDFT, L2OPTIMALMODEAMPLITUDES

% Copyright 2015 under BSD license (see LICENSE file).

  import koopman.*

  % We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
  % column-by-column
  [InputSnapshots,OutputSnapshots] = debias(Snapshots, varargin{:});

  [Q,R] = qr(InputSnapshots,0);

  S = pinv(R) * Q' * OutputSnapshots;

  [X, lambdas] = eigs(S, size(S,1)-2);
  lambdas = diag(lambdas);

  %% Calculate modes
  Modes = InputSnapshots * X;
  lambdas = log(lambdas);
  if nargin >= 2 && ~isempty(dt)
    lambdas = lambdas/dt;
  end

  %%
  % Sort modes according to their optimal L2 contributions
  [~,lambdas, Modes, Amps] = sortmodes( lambdas, Modes, Snapshots, dt );
