function [lambdas, Modes, Amps] = DMD(Snapshots, dt, db)
%DMD Compute Koopman modes by "exact" Dynamic Mode Decomposition of Tu et al.
%
% This is the algorithm by Tu, Jonathan H., Clarence W. Rowley, Dirk
% M. Luchtenburg, Steven L. Brunton, and J. Nathan Kutz. 2013. “On Dynamic
% Mode Decomposition: Theory and Applications.” Journa of Computational
% Dynamics, doi:10.3934/jcd.2014.1.391.
%
% [lambdas, Modes, Amps] = DMD( Snapshots, dt )
%    Compute DMD of data in Snapshots matrix. Columns of Snapshots are
%    measurements taken dt apart.
%
%    lambdas -- list of complex Dynamic Mode frequencies, real part is the
%    decay rate, imaginary part (angular) frequency.
%    Modes  -- each L2-normalized column of the matrix is a Dynamic Mode, corresponding
%    to the lambdas at the same index.
%    Amps   -- optimal L2 amplitudes used to sort the modes in descending order
%
% [lambdas, Modes, Amps] = DMD( ..., db ) If set to true, debias first
%    "de-biases" the data using truncation of SVD modes of the snapshot
%    matrix, according to:
%
%    Hemati, Maziar S., and Clarence W. Rowley. 2015. “De-Biasing the Dynamic
%    Mode Decomposition for Applied Koopman Spectral Analysis.”
%    arXiv:1502.03854 [physics], February. http://arxiv.org/abs/1502.03854.
%
%    By default, debias uses all non-negligible SVD directions
%    directions. Alternatively, user can request a specific number of
%    directions by passing an integer through debias.
%
%
% See also DMD_DUKE, KDFT, L2OPTIMALMODEAMPLITUDES

% Copyright 2015 under BSD license (see LICENSE file).

  import koopman.*

  % We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
  % column-by-column
  [InputSnapshots,OutputSnapshots] = debias(Snapshots, db);

  %% "Exact" DMD
  [Ux, Sx, Vx] = svd( InputSnapshots, 'econ' );

  % invert the Sigma matrix, skipping over almost-zero elements
  SxInv = pinv(Sx);

  %% Eigenvectors of Atilde will give Koopman modes
  Atilde = Ux' * OutputSnapshots * Vx * SxInv;

  [w, lambdas] = eigs(Atilde, size(Atilde,1)-2);
  lambdas = diag(lambdas);

  %% Calculate modes
  %  Modes = Snapshots(:,2:end) * Vx * SxInv * w * diag(1./lambdas);
  Modes = Ux * w;

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
