function [lambda, Modes] = DMD_Duke(Snapshots, dt, varargin)
%DMD_DUKE Compute Koopman modes by Dynamic Mode Decomposition algorithm by Duke et al.
%
% This is the algorithm by Duke, Daniel, Julio Soria, and Damon
% Honnery. 2012. “An Error Analysis of the Dynamic Mode Decomposition.”
% Experiments in Fluids 52 (2): 529–42. doi:10.1007/s00348-011-1235-7.
%
% [lambda, Modes] = DMD_Duke( Snapshots, dt )
%    Compute DMD of data in Snapshots matrix. Columns of Snapshots are
%    measurements taken dt apart.
%
%    lambda -- list of complex Dynamic Mode frequencies, real part is the
%    decay rate, imaginary part (angular) frequency.
%    Modes  -- each column of the matrix is a Dynamic Mode, corresponding
%    to the lambda at the same index.
%
%    lambda and Modes are sorted by l2-norm of columns of Modes, in
%    descending order.
%
% [lambda, Modes] = DMD_Duke( ..., db ) If db set to true, first
%    "de-bias" using Hemati, Rowley procedure. Alternatively, if db is a
%    positive integer, use db modes to debias.
%
%
% See also DMD, KDFT

% Copyright 2015 under BSD license (see LICENSE file).

  % We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
  % column-by-column
  [InputSnapshots,OutputSnapshots] = debias(Snapshots, varargin{:});

  [Q,R] = qr(InputSnapshots);

  S = pinv(R) * Q' * OutputSnapshots;

  [X, lambda] = eigs(S, size(S,1)-2);
  lambda = diag(lambda);

  %% Calculate modes
  Modes = InputSnapshots * X;
  lambda = log(lambda);
  if nargin >= 2 && ~isempty(dt)
    lambda = lambda/dt;
  end

  [lambda, Modes] = sortmodes( lambda, Modes );
  Modes = bsxfun( @rdivide, Modes, koopman.columnNorm(Modes) );
