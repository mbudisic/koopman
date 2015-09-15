function DukeCalibration(isNoisy)
%DUKECALIBRATION Compute Koopman modes using several techniques on
%synthetic data.
%
% This function implements the data set similar to
% Duke, Daniel, Julio Soria, and Damon Honnery. 2012. “An Error Analysis of
% the Dynamic Mode Decomposition.” Experiments in Fluids 52 (2):
% 529–42. doi:10.1007/s00348-011-1235-7.

% DUKECALIBRATION(ISNOISY)
% Data set contains an exponential spatial and temporal shape, with added
% noise (SNR=100) (unless FALSE) is passed as the argument
%
% It then computes Koopman modes using exact, Duke, and DFT algorithms and
% plots them.
%
% This function should be interpreted as a "sanity" check for Koopman mode
% techniques.
%

% Copyright 2015 under BSD license (see LICENSE file).

import koopman.*

%%
% Set up the domains
xPoints = 201; % points on space domain
tPoints = 501; % points on time domain

% Generate Duke Synthetic Data set
[U, t, x] = DukeSynthetic(xPoints, tPoints, 20i, 1+5i);

dt = t(2)-t(1); % sample time
dx = x(2)-x(1); % sample step

% Add multiplicative noise
if nargin == 1 && ~isNoisy
  disp('Noiseless')
else
  disp('Adding noise')
  NSR = 5/100; % noise to signal ratio
  Noise = (2*rand(size(U)) - 1) * NSR;
  U = U .* (1 + Noise);
end


figure
pcolor( t, x, U );
xlabel('Time t');
ylabel('Space x');
shading interp
title('Synthetic Duke data')
colorbar

figure
subplot(1,2,1)
tPeak = getspectrum( U(end,:), dt );
xPeak = getspectrum( U(:,end), dx );

getspectrum( U(end,:), dt );

fprintf('FFT-detected time ang. freq: %f\n', tPeak(1,1));
fprintf('FFT-detected space ang. freq: %f\n', xPeak(1,1));

Shape = U(:,end);
Shape = Shape/max(abs(Shape));
Nmd = 5;

[U, Mean] = removemean(U);

tic
[lambda_u1, Phi_u1, Amp_u1] = DMD( U, dt, 20 );
toc
tic
[lambda_u2, Phi_u2, Amp_u2] = DMD_Duke( U, dt, 20 );
toc
tic
[lambda_u3, Phi_u3, Amp_u3] = KDFT( U, dt  );
toc

normalize = @(v)v/max(abs(v));



subplot(1,2,2);
x = x.';

h = plot(x,U(:,1),'LineWidth',3 );
h.DisplayName = 'Data';
hold all;

plotMode( x, Amp_u1(1)*Phi_u1(:,1), 'Exact DMD' );
plotMode( x, Amp_u2(1)*Phi_u2(:,1), 'Duke DMD' );
plotMode( x, Amp_u3(1)*Phi_u3(:,1), 'KDFT' );

legend('Location','Best');

disp('Exact DMD:')
Amp_u1(1:Nmd).'
lambda_u1(1:Nmd)

disp('Duke DMD:')
Amp_u2(1:Nmd).'
lambda_u2(1:Nmd)

disp('KDFT:')
Amp_u3(1:Nmd).'
lambda_u3(1:Nmd)

end

function h = plotMode( x, z, name )
%PLOTMODE Plot the complex mode.

  validateattributes(x, {'numeric'},{'column','finite','nonnan'});
  validateattributes(z, {'numeric'},{'column','finite',...
                      'nonnan','numel',numel(x)});

  % multiply by 2 to compensate for the conjugate mode
  h = plot( x, 2*real(z) );
  h.DisplayName = name;

end
