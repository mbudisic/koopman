function AMP = L2optimalModeAmplitudes( lambdas, Modes, Snapshots, t )
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
%

import koopman.*

% validate inputs
validateattributes( lambdas, {'numeric'}, {'column'} )
validateattributes( t, {'numeric'}, {'row'} )
validateattributes( Snapshots, {'numeric'},...
                    {'2d','ncols',numel(t)});
validateattributes( Modes, {'numeric'},...
                    {'2d','nrows',size(Snapshots,1)});

% scale modes to unit norm
%Modes = bsxfun( @rdivide, Modes, columnNorm(Modes) );

% generalized Vandermonde matrix
Vand = exp( bsxfun( @times, lambdas, t ) );

[U,Sigma,V] = svd(Snapshots);

Y = U'*Modes;
P = (Y'*Y) .* conj(Vand*Vand');
q = conj( diag(Vand*V*Sigma'*Y) );

Amp = pinv(P)*q;

% Amp = abs(Alpha);
% Angles = Alpha ./ Amp;
% PhasedModes = bsxfun( @times, Modes, conj(Angles(:).') );
