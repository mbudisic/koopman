function CompFun = DemoKoopmanModes(TCF, noNoise)
%%DEMOKOOPMANMODES Demonstrate Koopman mode calculation.
%
% This function omputes Koopman modes using several different techniques
% (exact DMD, Duke DMD, stabilized snapshot DMD, and Koopman DFT).
% from a synthetic data set, similar to the one used in
% Duke, Daniel, Julio Soria, and Damon Honnery. 2012. “An Error Analysis of
% the Dynamic Mode Decomposition.” Experiments in Fluids 52 (2):
% 529–42. doi:10.1007/s00348-011-1235-7.
% The function can be interpreted as a "sanity check" for implementations
% of Koopman mode decompositions.
%
% The data set is a complex exponential in one spatial and temporal
% dimension, with angular time-frequency 20, and angular space-frequency 5,
% and growth rate 1.
%
% DEMOKOOPMANMODES Compute Koopman modes for a data set containing an
% exponential spatial and temporal shape, with added multplicative noise
% (Signal-Noise-Ratio = 20).
%
% DEMOKOOPMANMODES(TCF) Compute Koopman modes for a data set containing an
% exponential spatial and temporal shape, with explicit complex temporal
% frequency given.
%
% DEMOKOOPMANMODES(TCF, NONOISE) Compute Koopman modes for a data set containing an
% exponential spatial and temporal shape, with explicit complex temporal
% frequency given. Do not add noise if NONOISE = true. Default TCF is set
% by TCF=[]

% See DUKESYNTHETIC

% Copyright 2015 under BSD license (see LICENSE file).

import koopman.*

% plot Nmd most dominant modes
Nmd = 10;

%% Generate Duke Synthetic Data set
if nargin < 1 || isempty(TCF)
  TCF = -0.1 + 21i;
else
  validateattributes(TCF,{'numeric'},{'scalar','finite','nonnan'})
end

fprintf('Complex time frequency: %.1f  + i %.1f\n', real(TCF), imag(TCF) );

[U, t, x] = DukeSynthetic('TimeComplexFrequency', TCF, ...
                          'SpaceComplexFrequency', 1+5i);
% compute time and space step sizes
dt = t(2)-t(1);
dx = x(2)-x(1);

%%
% If requested, add multiplicative noise
if nargin > 1 && noNoise
  disp('Noiseless')
else
  disp('Adding noise')
  NSR = 10/100; % noise to signal ratio
  Noise = (2*rand(size(U)) - 1) * NSR;
  U = U .* (1 + Noise);
end

%%
% Plot the time-space false color plot of data
subplot(2,1,1);
pcolor( t, x, U );
xlabel('Time t');
ylabel('Space x');
shading interp
title('Synthetic Duke data set')
c = colorbar('South');
c.Position(4) = c.Position(4)/4;

%%
% Compute space and time power spectra for the last snapshot in the data
% set and label peaks in the spectrum
subplot(2,2,3)
getPowerSpectrum( U(end,:), dt );
title('Time FFT')

subplot(2,2,4)
getPowerSpectrum( U(:,end), dx );
title('Space FFT')

%% Compute Koopman Modes

%%
% Remove mean from data
[U, Mean] = removemean(U);
fprintf('Removed data mean (ranged in interval [%f, %f])\n', ...
        min(Mean), max(Mean) );

%
% CompFun stores calls of functions to compute Koopman modes
CompFun = struct([]);

CompFun(end+1).Eval = @(Data)DMD( Data, dt, 20 );
CompFun(end).Name = 'Exact DMD (de-biased)';

CompFun(end+1).Eval = @(Data)DMD_Duke( Data, dt, 20 );
CompFun(end).Name = 'Duke DMD (de-biased)';

CompFun(end+1).Eval = @(Data)DMD_Snapshot( Data, dt, 20 );
CompFun(end).Name = 'Snapshot DMD (de-biased)';

CompFun(end+1).Eval = @(Data)KDFT( Data, dt);
CompFun(end).Name = 'Koopman DFT';

%%
% Compute modes
for k = 1:numel(CompFun)
  fprintf('Calculating %s.\n',CompFun(k).Name);
  tic
  [CompFun(k).Spectrum, ...
   CompFun(k).Modes, ...
   CompFun(k).Amplitudes] = CompFun(k).Eval(U);
  toc
end

%% Plotting

figure('Name', 'Results');

%%
% Plot mode shapes
subplot(2,2,1);

% Plot data
x = x(:);
step = numel(t);

truth = U(:,step);
h = plot(x,truth,'k-','LineWidth',2 );
h.DisplayName = 'Data';
hold all;
axis manual; % fix axis according to data

% Plot modes
for k = 1:numel(CompFun)
  fprintf('Plotting %s.\n',CompFun(k).Name);

  % multiply by 2 to compensate for the conjugate mode
  y = 2*real(...
      CompFun(k).Amplitudes(1)*...
      CompFun(k).Modes(:,1)*...
      exp(CompFun(k).Spectrum(1)*t(step))...
      );
  d = y - truth;

  subplot(2,2,1);
  h = plot( x, y, '.-', 'LineWidth',1);
  h.DisplayName = CompFun(k).Name;

  subplot(2,2,3);
  h = plot( x, d, '.-', 'LineWidth',1);
  h.DisplayName = CompFun(k).Name;
  hold all
end

subplot(2,2,1);
title({'Mode shapes compared to data',...
       sprintf('in time step %d/%d', step, numel(t))});
hl = legend('Location','Best');
hl.FontSize = 7;

subplot(2,2,3);
title({'Difference betweend data and mode shape',...
       sprintf('in time step %d/%d', step, numel(t))});
hl = legend('Location','Best');
hl.FontSize = 7;

hold off;

%%
% Plot spectra
subplot(1,2,2);
h = gobjects(numel(CompFun),1);
markers = {'o','^','v','s','d'};
for k = 1:numel(CompFun)
  MySpectrum = CompFun(k).Spectrum(1:Nmd);
  MyAmps = abs(CompFun(k).Amplitudes(1:Nmd));

  Xs = real(MySpectrum);
  Ys = imag(MySpectrum);
  Ss = MyAmps*25/max(MyAmps);

  h(k) = scatter( Xs, Ys, Ss );

  h(k).DisplayName = CompFun(k).Name;

  %fix Matlab bug that assigns last color to all previously plotted scatter
  %points
  h(k).MarkerFaceColor = h(k).CData(1,:);
  h(k).MarkerEdgeColor = 'none';
  h(k).Marker = markers{k};
  hold all;
end
axis( [ [-2,2]*max([abs(real(TCF)),.1]),...
        [-2,2]*abs(imag(TCF)) ] )
hold off;
hl = legend('Location','Best');
hl.FontSize = 7;
xlabel('Decay Rate');
ylabel('Frequency');
title({sprintf('Dominant (%d) Koopman evalues',Nmd),...
       '(size is mode amplitude)'});

end
