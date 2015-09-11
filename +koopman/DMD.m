function [lambda, Modes] = DMD(Snapshots, dt, Nmd, db)
%DMD Compute Koopman modes by "exact" Dynamic Mode Decomposition of Tu et al.
%
% This is the algorithm by Tu, Jonathan H., Clarence W. Rowley, Dirk
% M. Luchtenburg, Steven L. Brunton, and J. Nathan Kutz. 2013. “On Dynamic
% Mode Decomposition: Theory and Applications.” Journa of Computational
% Dynamics, doi:10.3934/jcd.2014.1.391.
%
% [lambda, Modes] = DMD( Snapshots, dt, Nmd )
%    Compute DMD of data in Snapshots matrix. Columns of Snapshots are
%    measurements taken dt apart. Nmd is the number of modes to be
%    computed.
%
%    lambda -- list of complex Dynamic Mode frequencies, real part is the
%    decay rate, imaginary part (angular) frequency.
%    Modes  -- each column of the matrix is a Dynamic Mode, corresponding
%    to the lambda at the same index.
%
%    lambda and Modes are sorted by l2-norm of columns of Modes, in
%    descending order.
%
% [lambda, Modes] = DMD( ..., db ) If set to true, debias first
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

% Copyright 2015 under BSD license (see LICENSE file).

  % We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
  % column-by-column
  [InputSnapshots,OutputSnapshots] = debias(Snapshots, db);

  %% "Exact" DMD
  [Ux, Sx, Vx] = svd( InputSnapshots, 'econ' );

  % invert the Sigma matrix, skipping over almost-zero elements
  SxInv = pinv(Sx);

  %% Eigenvectors of Atilde will give Koopman modes
  Atilde = Ux' * OutputSnapshots * Vx * SxInv;

  [w, lambda] = eigs(Atilde, [], Nmd, 'LR');
  lambda = diag(lambda);

  %% Calculate modes
  %  Modes = Snapshots(:,2:end) * Vx * SxInv * w * diag(1./lambda);
  Modes = Ux * w;

  lambda = log(lambda);
  if ~isempty(dt)
    lambda = lambda/dt;
  end

  [lambda, Modes] = sortmodes( lambda, Modes );