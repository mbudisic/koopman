function Modes = normalize(Modes)
%NORMALIZE Normalize columns of Modes matrix by their L2 norms
%
% MODES = NORMALIZE(MODES)
%

NormOfModes = sqrt( sum( abs(Modes).^2, 1 ) / size(Modes,1) );
Modes = bsxfun( @rdivide, Modes, NormOfModes );
