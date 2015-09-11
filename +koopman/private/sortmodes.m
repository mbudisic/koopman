function [lambda, Modes] = sortmodes( lambda, Modes, sorttype );
%SORTMODES Sort modes and frequencies according to the argument sorttype.
%
%  Currently, sorttype is ignored and L2 vector norm of modes is used for
%  sorting.
%

% Copyright 2015 under BSD license (see LICENSE file).

% assert the sizes match
assert( numel(lambda) == size(Modes,2), ...
        ['Number of frequencies and modes has to match']);

% compute the 2-norm of Modes as vectors
NormOfModes = sqrt( sum( abs(Modes).^2, 1 ) / size(Modes,1) );

% sort the Norms in descending order, and use the order to sort frequencies
% and modes themselves
[NormOfModes,Order] = sort( NormOfModes, 2, 'descend' );
lambda = lambda(Order);
Modes = Modes(:,Order);

disp('ModeNorms')
NormOfModes(1:min(5, numel(lambda)))
