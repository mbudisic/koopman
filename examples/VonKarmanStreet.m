%% KOOPMAN DECOMPOSITION OF VON KARMAN STREET
% This is a demo of Koopman decomposition of a fluid flow data.
%%

% Copyright 2017 under BSD license (see LICENSE file).

function VonKarmanStreet( recompute )

vonkarmanurl='https://www.dropbox.com/s/p0cl8t7q9l2qwe4/VonKarmanStreet.mat?dl=1';
vonkarmanlocal='VonKarmanStreet.mat';

%% Load the velocity field
if ~exist(vonkarmanlocal,'file')
  warning('%s file missing.', vonkarmanlocal);
  fprintf(['Press <Enter> to download the file (ca 170MB) from Dropbox to this folder.\n URL: %s \n Otherwise press ' ...
           'Ctrl-C to stop.\n'], vonkarmanurl);
  pause
  disp('Downloading (this may take a few minutes, depending on your connection)');
  websave(vonkarmanlocal, vonkarmanurl );
end

assert( exist(vonkarmanlocal,'file') == 2, '%s does not exist', vonkarmanlocal );

fprintf('Loading %s',vonkarmanlocal);
load(vonkarmanlocal);
disp(DESCRIPTION)
disp(README)



%% Compute speed and vorticity
% These quantities will be used as measurements
disp('Computing speed and vorticity (observables)...');
%
% Evaluate discrete derivatives used to compute vorticity
Ux = diff(U, 1, 3 );
Uy = diff(U, 1, 2 );

% Pad to preserve data size
Ux(:,:,end+1,:) = Ux(:,:,end,:);
Uy(:,end+1,:,:) = Uy(:,end,:,:);

% vorticity
Vorticity = Ux(2,:,:,:) - Uy(1,:,:,:);

% speed
Speed = sqrt(sum(U.^2,1));

% stack data
Data = cat(1, Speed, Vorticity );
[nVars, nRows, nCols, nSteps] = size(Data);
% Data is a 4-d array

% plot the first step
plotData(Data, 1)



%% Format data to suit Koopman toolbox
disp('Formatting to suit Koopman toolbox...');

% Each column in Snapshots contains all the data recorded at a single timestep
Snapshots = reshape(Data, nVars*nRows*nCols, nSteps );
% Snapshots is now a 2d array (matrix)

%%
% Remove the mean
Snapshots = koopman.removemean(Snapshots);

%% Compute Koopman modes

% Use persistent variables to avoid recomputing Modes,
% e.g., when just re-plotting data
persistent Spectrum
persistent Modes
persistent Amps

% If no modes have been computed, or recomputation is explictly requested,
% compute the Koopman modes
if isempty(Spectrum) || (nargin >= 1 && recompute)

  %%
  % Main invocation of the computation
  [Spectrum, Modes, Amps] = koopman.DMD( Snapshots, 1, true);

  %%

  % Reshape modes back into Data-style 4-d format
  nModes = numel(Spectrum);
  Modes = real(reshape(Modes, [nVars,nRows,nCols,nModes]));
end

%% Plot the spectrum
figure('Name','Spectrum')
koopman.plotSpectrum( Spectrum, Amps )

%%
% The size of the marker indicates the L2-optimal magnitude of the mode, and
% its color indicates the L2-optimal phase of the mode

%% Plot the modes
figure('Name','Modes');
plotModes(4, Spectrum, Amps, Modes)

%% Auxiliary plotting functions

%%% PLOTDATA
% Plot one step in data matrix
function plotData( Data, step )

figure('Name','Data');
subplot(2,1,1);
pcolor(squeeze(Data(1,:,:,step)))
shading interp
xlabel('X');
ylabel('Y');
title(sprintf('Speed at step %d',step) )
axis equal tight
set(gca,'XTick',[],'YTick',[])

subplot(2,1,2);
pcolor(squeeze(Data(2,:,:,step)))
shading interp
xlabel('X');
ylabel('Y');
title(sprintf('Vorticity at step %d',step) )
axis equal tight
set(gca,'XTick',[],'YTick',[])

%%% PLOTMODES
% Plot speed and vorticity modes in columns

function plotModes( nPlotted, Spectrum, Amps, Modes )

%%
% Visualize only modes with positive angles -- to avoid replication between
% conjugate pairs of modes

sel = angle(Spectrum) >= 0;
Spectrum = Spectrum(sel);
Amps = Amps(sel);
Modes = Modes(:,:,:,sel);

%%
% Split up modes for speed and vorticity
SpeedModes = squeeze(Modes(1,:,:,:));
VortModes = squeeze(Modes(2,:,:,:));

for k = 1:nPlotted

  %%
  % Plot Speed mode
  subplot(nPlotted, 2, 2*k-1);
  pcolor( SpeedModes(:,:,k) ); shading interp;
  title({sprintf('%d: Growth=%.1f, Freq=%.1f', k, real(Spectrum(k)),...
                 imag(Spectrum(k)) ), ...
         'Speed'});
  set(gca,'XTick',[],'YTick',[])
  axis equal tight

  %%
  % Plot Vorticity mode
  subplot(nPlotted, 2, 2*k);
  pcolor( VortModes(:,:,k) ); shading interp;

  title({sprintf('log_{10}(Amplitude)=%.1f, Phase=%.1f \\pi', ...
                 log10(abs(Amps(k))),...
                 angle(Amps(k))/k ), ...
         'Vorticity'});
  set(gca,'XTick',[],'YTick',[])
  axis equal tight
end
