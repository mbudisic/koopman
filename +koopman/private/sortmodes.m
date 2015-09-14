function varargout = sortmodes( lambdas, Modes, Snapshots, T )
%SORTMODES Sort modes and frequencies based on their contribution to snapshots.
%
%  IND = SORTMODES( LAMBDAS, MODES, SNAPSHOTS, T )
%  Compute indices used to sort modes. Modes are L2-normalized and sorted
%  in decreasing order according to their L2-contribution in reconstruction
%  of SNAPSHOTS. T contains times at which SNAPSHOTS columns were
%  obtained. If T is a scalar, uniform sampling at time step T is assumed.
%
%  [IND, LAMBDAS, MODES, AMPS] = SORTMODES( ... )
%  As before, but also return sorted frequencies, normalized modes, amplitudes,
%

% Copyright 2015 under BSD license (see LICENSE file).

%%
% Validate inputs
assert( nargin == 4, 'Specify all 4 arguments');

validateattributes(Modes, {'numeric'}, {'2d','finite'});
D = size(Modes, 1);
M = size(Modes, 2);

validateattributes(Snapshots, {'numeric'}, {'2d','finite','nrows',D});
N = size(Snapshots,2);

validateattributes(lambdas, {'numeric'}, {'column','finite','numel',M});

if isscalar(T)
  validateattributes(T, {'numeric'}, {'finite','positive','scalar'});
  T = T*(0:(N-1));
else
  validateattributes(T, {'numeric'}, {'increasing','finite','row'});
end

%%
% Normalize columns of Modes matrix by their L2 norms
NormOfModes = sqrt( sum( abs(Modes).^2, 1 ) / size(Modes,1) );
Modes = bsxfun( @rdivide, Modes, NormOfModes );

%%
% Compute L2 optimal amplitudes
Amps = koopman.L2optimalModeAmplitudes( lambdas, Modes, Snapshots, T );
validateattributes(Amps,{'numeric'},{'column', 'numel', M})

%%
% sort the Norms in descending order, and use the order to sort frequencies
% and modes themselves
[Amps,idx] = sort( Amps, 1, 'descend' );
lambdas = lambdas(idx);
Modes = Modes(:,idx);

%%
% Assign outputs
if nargout > 0
  varargout = {};

  if nargout >= 1
    varargout{1} = idx;
  end
  if nargout >= 2
    varargout{2} = lambdas;
  end
  if nargout >= 3
    varargout{3} = Modes;
  end
  if nargout >= 4
    varargout{4} = idx;
  end
end

end
