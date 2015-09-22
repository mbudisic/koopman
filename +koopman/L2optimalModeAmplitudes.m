function Amp = L2optimalModeAmplitudes( lambdas, Modes, Snapshots, t )
%L2OPTIMALMODEAMPLITUDES Compute optimal L2 reconstruction of data by modes.
%
% Use L2-residual minimization to compute optimal (complex) amplitudes that
% reconstruct SNAPSHOTS as a sum of MODES oscillating at frequencies LAMBDA.
%
% See in particular Sec IIB in Jovanović, Mihailo R., Peter J. Schmid, and
% Joseph W. Nichols. 2014. “Sparsity-Promoting Dynamic Mode Decomposition.”
% Physics of Fluids (1994-Present) 26 (2): 024103.
%
% AMP = L2OPTIMALMODEAMPLITUDES( LAMBDAS, MODES, SNAPSHOTS, T )
%   Returns a column-vector of amplitudes.
%
% See also DMD_DUKE, KDFT, DMD, DMD_SNAPSHOT

% Copyright 2015 under BSD license (see LICENSE file).


import koopman.*

% validate inputs
validateattributes( lambdas, {'numeric'}, {'column'} )
validateattributes( t, {'numeric'}, {'row'} )
validateattributes( Snapshots, {'numeric'},...
                    {'2d','ncols',numel(t)});
validateattributes( Modes, {'numeric'},...
                    {'2d','nrows',size(Snapshots,1)});

% generalized Vandermonde matrix
Vand = exp( bsxfun( @times, lambdas, t ) );

%%
% Solve L2 optimization, see Jovanovic et al (2014)
[U,Sigma,V] = svd(Snapshots,'econ');

Y = U'*Modes;
P = (Y'*Y) .* conj(Vand*Vand');
q = conj( diag(Vand*V*Sigma'*Y) );

Amp = pinv(P)*q;
