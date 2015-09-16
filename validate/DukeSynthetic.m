function [U, t, x, peakTimeFFT, peakSpaceFFT] = DukeSynthetic( varargin )
%%DUKESYNTHETIC Generate a synthetic 1-d linear instability field
% from D. Duke, J. Soria, and D. Honnery, “An error analysis of the dynamic
% mode decomposition,” Exp. Fluids 52, 529–542 (2012).
%
% [U, t, x] = DUKESYNTHETIC( NAME, VALUE, ... )
% Evaluate the evolution of an exponential space-time mode
% exp( SpaceFrequency x ) * ( TimeFrequency t )
% where both Space and Time frequencies are complex numbers, whose real
% parts are growth/decay rates, and imaginary parts oscillation
% frequencies.
%
% Function takes the following optional named parameters (default value in
% parentheses):
%
%   TimeLength -- length of time interval (default: 5)
%   SpaceWidth -- width of the space domain (default: 3)
%   TimePoints -- number of points on the time domain (default: 501)
%   SpacePoints -- number of points in the space domain (default: 201)
%   TimeComplexFrequency -- complex time frequency (default: 20i)
%   SpaceComplexFrequency -- complex space frequency (default: 0.5 + 5i)
%
% The outputs are:
% U -- data matrix, each row corresponds to a time-slice, each column to a
%      time trace of the mode evolution
% t -- time axis
% x -- space axis
% peakTimeFFT -- expected peak of the time FFT
% peakSpaceFFT -- expected peak of the space FFT
%
% If no output arguments are requested, a false-color plot of the data is
% produced.
%
% See also DEMOKOOPMANMODES

% Copyright 2015 under BSD license (see LICENSE file).

% parse arguments
parser = inputParser;

%%
% validators for parameter types

positiveScalar = @(n)isscalar(n) && ...
               isnumeric(n) && ...
               isfinite(n) && ...
               n > 0;

complexScalar = @(n)isscalar(n) && ...
               isnumeric(n) && ...
               isfinite(n);

parser.addParameter('TimeLength', 5, positiveScalar );
parser.addParameter('SpaceWidth', 3, positiveScalar );
parser.addParameter('TimePoints', 501, positiveScalar );
parser.addParameter('SpacePoints', 201, positiveScalar );
parser.addParameter('TimeComplexFrequency', 20i, complexScalar);
parser.addParameter('SpaceComplexFrequency', 0.5+5i, complexScalar);
parser.parse(varargin{:});

p = parser.Results;

%%
% Extract real and imaginary parts
Omega = imag(p.TimeComplexFrequency);
Sigma = real(p.TimeComplexFrequency);

Kappa = imag(p.SpaceComplexFrequency);
Gamma = real(p.SpaceComplexFrequency);

% Mode function
Mode = @(t,x)real( exp( p.TimeComplexFrequency*t ) .* ...
                   exp( p.SpaceComplexFrequency*x ) );

%%
% Compute theoretical peaks of the power spectrum
ZetaT = sqrt( Sigma^2/(Omega^2 + Sigma^2) );
ZetaX = sqrt( Gamma^2/(Kappa^2 + Gamma^2) );
peakTimeFFT = Omega/sqrt(1-ZetaT.^2);
peakSpaceFFT = Kappa/sqrt(1-ZetaX.^2);

fprintf('Time Peak: %f\n',peakTimeFFT)
fprintf('Space Peak: %f\n',peakSpaceFFT)

%%
% Evaluate the mode on the time-space grid
x = linspace(0,p.SpaceWidth,p.SpacePoints);
t = linspace(0,p.TimeLength,p.TimePoints);
U = bsxfun( Mode, t, x' );

if nargout < 1
  pcolor( t, x, U );
  xlabel('Time t');
  ylabel('Space x');
  shading interp
end
