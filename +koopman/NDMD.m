function [lambda, M] = NDMD( snapshots, t, Nmd )
%NDMD Compute DMD for nonuniformly sampled data.
%
% Uses algorithm by Gueniat et al 2015
%
% NDMD( snapshots, t, Nmd )
%
% snapshots -- data matrix, each column is a time slice of data
% t -- vector of times at which data was recorded
% Nmd -- number of modes to retain

  validateattributes(t, {'numeric'},{'vector','increasing'});
  tRef = t(1);
  t = t-tRef;
  N = numel(t);
  t = t(:)'; % make into a row-vector

  if ~isscalar(Nmd)
    guess = Nmd;
    Nmd = numel(guess);
  else
    guess = complex( rand([Nmd,1]),rand([Nmd,1]) );
  end

  tM = repmat( t, [Nmd, 1] );

  dt = min(diff(t));

  validateattributes( snapshots, {'numeric'}, {'2d','ncols', numel(t)} );

  np = size(snapshots, 1); % number of observables

  [~,R] = qr( snapshots, 0 );

  function L = Lambda( z )
    z = repmat( exp(z(:)), [1,N] );
    L = z .^ tM;
    L = [L; conj(L)];
  end

  function Res = matrixRes( z )
    Res = R*( eye(N) - peye( Lambda(z) ) );
  end

  normRes = @(z) norm( matrixRes(z), 'fro').^2/2;

  disp('Guess for NDMD:')
  lambda = ComplexNLSQ( @matrixRes, guess )
  lambda = ComplexOptimize( normRes, guess )
  normRes(lambda)

  M = snapshots*pinv(Lambda(lambda));

end
