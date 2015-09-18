function [lambdas, Modes, Amps] = DMD_Snapshot(Snapshots, dt, varargin)
%DMD Compute Koopman modes by "exact" Dynamic Mode Decomposition of Tu et al.
%
% This is the stabilized snapshots algorithm by Chen, Kevin K, Jonathan H Tu,
% and Clarence W Rowley. 2012. “Variants of Dynamic Mode Decomposition:
% Boundary Condition, Koopman, and Fourier Analyses.” Journal of Nonlinear
% Science, April. doi:10.1007/s00332-012-9130-9.
%
% The algorithm uses the SVD decomposition of the square of the input snapshot
% matrix to construct the DMD matrix.
%
% [lambdas, Modes, Amps] = DMD_Snapshot( Snapshots, dt )
%    Compute DMD of data in Snapshots matrix. Columns of Snapshots are
%    measurements taken dt apart.
%
%    lambdas -- list of complex Dynamic Mode frequencies, real part is the
%    decay rate, imaginary part (angular) frequency.
%    Modes  -- each L2-normalized column of the matrix is a Dynamic Mode, corresponding
%    to the lambdas at the same index.
%    Amps   -- optimal L2 amplitudes used to sort the modes in descending order
%
% [lambdas, Modes, Amps] = DMD_Snapshot( ..., db ) If set to true, debias first
%    "de-biases" the data using truncation of SVD modes of the snapshot
%    matrix, according to:
%
%    Hemati, Maziar S., and Clarence W. Rowley. 2015. “De-Biasing the Dynamic
%    Mode Decomposition for Applied Koopman Spectral Analysis.”
%    arXiv:1502.03854 [physics], February. http://arxiv.org/abs/1502.03854.
%
%    By default, db = true, which uses all non-negligible SVD directions
%    directions. Alternatively, user can request a specific number of
%    directions by passing an integer through debias.
%
% The function returns both conjugate pairs of the modes and their frequencies.
%
% See also DMD, DMD_DUKE, KDFT, L2OPTIMALMODEAMPLITUDES

% Copyright 2015 under BSD license (see LICENSE file).

  import koopman.*

  % We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
  % column-by-column
  [InputSnapshots,OutputSnapshots] = debias(Snapshots, varargin{:});

  % Compute the modulus (square) of the input snapshot matrix
  % and its SVD
  Snapshots2 = InputSnapshots' * InputSnapshots;
  [W, Sigma2,~] = svd( Snapshots2, 'econ' );
  Sigma = sqrt(Sigma2);

  % Projection matrix
  U = InputSnapshots*W*pinv(Sigma);

  %%
  % Eigenvectors of Atilde will give Koopman modes
  Atilde = U' * OutputSnapshots * W * pinv(Sigma);
  [w, lambdas] = eigs(Atilde, size(Atilde,1)-2);
  lambdas = diag(lambdas);

  %% Calculate modes
  %  Modes = Snapshots(:,2:end) * Vx * SxInv * w * diag(1./lambdas);
  Modes = U * w;

  %%
  % Return complex arguments of lambdas, as they have physical
  % interpretations as DecayRate + 1j * Frequency
  lambdas = log(lambdas);
  if ~isempty(dt)
    lambdas = lambdas/dt;
  end

  %%
  % Sort modes according to their optimal L2 contributions
  [~,lambdas, Modes, Amps] = sortmodes( lambdas, Modes, Snapshots, dt );
