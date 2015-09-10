function DukeDataCalibrationTest

import koopman.*

xPoints = 201;
tPoints = 501;

[U1, t, x] = DukeSynthetic(xPoints, tPoints, -0.75 + 20i, 1+5i);
[U2, t, x] = DukeSynthetic(xPoints, tPoints, 0.75 + 15i, 1+10i);
U = U1;

NSR = 1/100; % noise to signal ratio
Noise = (2*rand(size(U)) - 1) * NSR;
%U = U .* (1 + Noise);


dt = t(2)-t(1);
dx = x(2)-x(1);

subplot(1,2,1);
pcolor( t, x, U );
xlabel('Time t');
ylabel('Space x');
shading interp

tPeak = getspectrum( U(end,:), dt );
xPeak = getspectrum( U(:,end), dx );

getspectrum( U(end,:), dt );

fprintf('FFT-detected time ang. freq: %f\n', tPeak(1,1));
fprintf('FFT-detected space ang. freq: %f\n', xPeak(1,1));

Shape = U(:,end);
Shape = Shape/max(abs(Shape));
Nmd = 4;

tic
[lambda_u1, Phi_u1] = DMD( U, dt, Nmd );
toc
tic
[lambda_u2, Phi_u2] = DMD_Duke( U, dt, Nmd );
toc
% tic
% [lambda_n, Phi_n] = NDMD( U, t, Nmd );
% toc

% normalize = @(v)v/max(abs(v));

% subplot(1,2,2);
% h = plot(x,[normalize(U(:,end)),...
%             real(normalize(Phi_u1(:,1))), ...
%             real(normalize(Phi_u2(:,2))),...
%             real(normalize(Phi_n(:,2)))] );
% h(1).DisplayName = 'Data';
% h(2).DisplayName = 'Exact DMD';
% h(3).DisplayName = 'Duke DMD';
% h(4).DisplayName = 'N-DMD';
% legend('Location','Best');

lambda_u1
lambda_u2
