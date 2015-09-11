function DukeCalibration(isNoisy)
%DUKECALIBRATION Compute Koopman modes using several techniques on
%synthetic data.
%
% This function implements the data set similar to
% Duke, Daniel, Julio Soria, and Damon Honnery. 2012. “An Error Analysis of the Dynamic Mode Decomposition.” Experiments in Fluids 52 (2): 529–42. doi:10.1007/s00348-011-1235-7.
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
if nargin == 1 && isNoisy
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
Nmd = 4;

[U, Mean] = removemean(U);

tic
[lambda_u1, Phi_u1] = DMD( U, dt, Nmd, true );
toc
tic
[lambda_u2, Phi_u2] = DMD_Duke( U, dt, Nmd, true );
toc
tic
[lambda_u3, Phi_u3] = KDFT( U, dt  );
toc

% tic
% [lambda_n, Phi_n] = NDMD( U, t, Nmd );
% toc

normalize = @(v)v/max(abs(v));

subplot(1,2,2);
h = plot(x,[normalize(U(:,end)),...
            real(normalize(Phi_u1(:,1))), ...
            real(normalize(Phi_u2(:,1))),...
            real(normalize(Phi_u3(:,1)))] );
h(1).DisplayName = 'Data';
h(2).DisplayName = 'Exact DMD';
h(3).DisplayName = 'Duke DMD';
h(4).DisplayName = 'KDFT';
legend('Location','Best');

lambda_u1(1:Nmd)
lambda_u2(1:Nmd)
lambda_u3(1:Nmd)
