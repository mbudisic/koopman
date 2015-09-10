function [U, tAx, xAx] = DukeSynthetic( xPoints, tPoints, TimeCF, SpaceCF )
%%DUKESYNTHETIC Generate a synthetic 1-d linear instability field
% from D. Duke, J. Soria, and D. Honnery, “An error analysis of the dynamic
% mode decomposition,” Exp. Fluids 52, 529–542 (2012).
%
% DukeSynthetic( xPoints, tPoints )
%
% Each row corresponds to a time-slice, each column to a time trace
%
% If called without return arguments, plot is produced.

% Parameters used in Gueniat et al. 2015
%Omega = 20; % time frequency
%Kappa = 10;  % space frequency

Omega = imag(TimeCF);
Sigma = real(TimeCF);

Kappa = imag(SpaceCF);
Gamma = real(SpaceCF);

ZetaT = sqrt( Sigma^2/(Omega^2 + Sigma^2) );
ZetaX = sqrt( Gamma^2/(Kappa^2 + Gamma^2) );

fprintf('Time Peak: %f\n',Omega/sqrt(1-ZetaT.^2))
fprintf('Space Peak: %f\n',Kappa/sqrt(1-ZetaX.^2))

U0 = 1;     % amplitude



xAx = linspace(0,2,xPoints);
tAx = linspace(0,5,tPoints);

[t,x] = meshgrid(tAx,xAx);

U = U0*imag( exp( ( Sigma - 1j*Omega )*t ) .* ...
        exp( (Gamma + 1j*Kappa) * x ) );

if nargout < 1
  pcolor( t, x, U );
  xlabel('Time t');
  ylabel('Space x');
  shading interp
end
