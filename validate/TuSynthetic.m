function [U, t, x, y] = TuSynthetic( varargin )
%TuSynthetic Generate a synthetic 2d field.
%
% The synthetic field consists of two spatial modes, each represented by a
% Gaussian blob, whose amplitude oscillates according to a sine function,
% with random offset of roughly 3%.
%
% Tu, Jonathan H., Clarence W. Rowley, J. Nathan Kutz, and Jessica
% K. Shang. 2014. “Spectral Analysis of Fluid Flows Using Sub-Nyquist-Rate PIV
% Data.” Experiments in Fluids 55 (9): 1–13. doi:10.1007/s00348-014-1805-6.
%

% Copyright 2015 under BSD license (see LICENSE file).

%% Parse input arguments
parser = inputParser;
positiveScalar = @(n)isscalar(n) && ...
               isnumeric(n) && ...
               isfinite(n) && ...
               n > 0;

parser.addParameter('TimeStep', 0.05, positiveScalar );
parser.addParameter('TimePoints', 8001, positiveScalar );
parser.addParameter('SpacePoints', 100, positiveScalar );

parser.addParameter('TimeFrequency1', 1.30, positiveScalar);
parser.addParameter('TimeFrequency2', 8.48, positiveScalar);
parser.parse(varargin{:});

p = parser.Results;

%% Construct the grid, time, and noise vectors
x = linspace(-2,2,p.SpacePoints);
y = linspace(-2,2,p.SpacePoints);
t = (0:p.TimePoints-1)*p.TimeStep;
noise = 0.1*rand(1, p.TimePoints);

%% Spatial modes
v1 = @(x,y) 2*exp( -(x-0.50).^2/(2*0.6^2) - (y-0.50).^2/(2*0.2^2) );
v2 = @(x,y)   exp( -(x+0.25).^2/(2*0.6^2) - (y-0.35).^2/(2*1.2^2) );

V1 = bsxfun(v1, x, y'); V1 = V1(:);
V2 = bsxfun(v2, x, y'); V2 = V2(:);

%% Time evolution
Osc1 = sin(p.TimeFrequency1*t);
Osc2 = sin(p.TimeFrequency2*t);

%% Combine space, time, and random data
S = bsxfun( @times, V1, Osc1 ) + bsxfun( @times, V2, Osc2 );
U = bsxfun( @plus, S, noise );

%% Animate data if no return arguments were requested
if nargout < 1

  U = reshape(U, p.SpacePoints,p.SpacePoints, p.TimePoints );

  %% Set up first image using the third slice (arbitrary)
  [~,h] = contourf( x, y, U(:,:,3) );
  xlabel('Space x');
  ylabel('Space y');
  shading interp
  caxis([-1,1]*max(abs(U(:))));
  caxis manual;

  %% Plot
  for k = 1:p.TimePoints
    h.ZData = U(:,:,k);
    pause(.1);
    title(sprintf('t = %.1f', t(k) ) );
  end

  %% Remove outputs so that they are not returned
  clear U x y t
end
