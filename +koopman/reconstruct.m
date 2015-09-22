function Output = reconstruct( t, Modes, Lambdas, Amplitudes, Baseline )
%RECONSTRUCT Reconstruct data using Koopman mode information.
%
% OUTPUT = RECONSTRUCT( T, MODES, LAMBDAS, AMPLITUDES, BASELINE )
% Reconstruct the output by appropriately summing Koopman modes.
% Inputs:
%   T          - row-vector of times at which to evaluate evolution
%   Modes      - complex matrix, each column is a Koopman mode
%   Lambdas    - complex column-vector, real parts are decay rates,
%                imaginary parts are angular frequencies of K. modes
%   Amplitudes - complex vector of amplitudes of Koopman modes
%   Baseline   - added to OUTPUT after K. mode reconstruction, e.g., mean
%                that was removed before processing, or a baseline fluid
%                flow; real 2d matrix, with either a single column, or
%                number of columns matching the number of time steps
%
% OUTPUT        - real reconstructed solution, each column corresponds
%                 to a timestep
%
% See also DMD, KDFT, DMD_DUKE, DMD_SNAPSHOT

% Copyright 2015 under BSD license (see LICENSE file).

%% Validate inputs

% Modes is a 2d matrix, determining sizes of all following inputs
validateattributes( Modes, {'numeric'},{'finite','nonnan','2d'} );

% number of measurement variables
D = size(Modes, 1);

% number of modes used
N = size(Modes, 2);

% times are given as a row-vector
validateattributes( t, {'numeric'},...
                    {'finite','nonnan','row','real'});

% complex frequencies are a column-vector, match the number of modes
validateattributes( Lambdas, {'numeric'},...
                    {'finite','nonnan','column','numel',N });

% complex amplitudes are a vector, match the number of modes
validateattributes( Amplitudes, {'numeric'},...
                    {'finite','nonnan','vector','numel',N });

% number of measurement variables in Baseline data has to match
% the number in modes
if nargin == 5
  validateattributes( Baseline, {'numeric'},...
                      {'real','finite','nonnan','nrows',D });
end

% generalized Vandermonde matrix
Vand = exp( bsxfun( @times, Lambdas, t) );

% combine the modes into the output
Output = real( Modes * diag(Amplitudes) * Vand );

% add the baseline if provided
if nargin == 5
  Output = bsxfun(@plus, Output, Baseline);
end
