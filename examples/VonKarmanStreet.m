function VonKarmanStreet( recompute )
% # Koopman decomposition of Von Karman street

% # Load the velocity field
disp('Loading data.')
load VonKarmanStreet.mat

% ##
% # We will use speed and vorticity as observables.
disp('Computing Speed and Vorticity')

% compute discrete derivatives and pad the matrices to recover the correct
% dimensions
Uy = diff(U, 1, 2 );
Uy(:,end+1,:,:) = Uy(:,end,:,:);

Ux = diff(U, 1, 3 );
Ux(:,:,end+1,:) = Ux(:,:,end,:);

% vorticity
W = Ux(2,:,:,:) - Uy(1,:,:,:);

% speed
V = sqrt(sum(U.^2,1));

%%
% Concatenate data
Data = cat(1, V, W );

%%
% Convert to Snapshots format -- stack all data at a single time step into
% its own column
nVars = size(Data,1);
nRows = size(Data,2);
nCols = size(Data,3);
nSteps=  size(Data,4);

disp('Computing Snapshots matrix');
Snapshots = reshape(Data, nVars*nRows*nCols, nSteps );

disp('Remove the mean')
Snapshots = koopman.removemean(Snapshots);

%%
% Compute Koopman modes
persistent Spectrum
persistent Modes
persistent Amps

if isempty(Spectrum) || (nargin >= 1 && recompute)
  disp('Computing Koopman modes');
  [Spectrum, Modes, Amps] = koopman.DMD( Snapshots, 1, true);
  % Reshape modes back into Data-style format
  nModes = numel(Spectrum);
  Modes = real(reshape(Modes, [nVars,nRows,nCols,nModes]));
else
  disp('Reusing already computed Koopman modes')
end

%%
% Plot the first frame of data
disp('Plotting data');
figure('Name','Data');
subplot(2,1,1);
pcolor(squeeze(Data(1,:,:,1)))
shading interp
xlabel('X');
ylabel('Y');
title('Speed')
axis equal tight
set(gca,'XTick',[],'YTick',[])

subplot(2,1,2);
pcolor(squeeze(Data(2,:,:,1)))
shading interp
xlabel('X');
ylabel('Y');
title('Vorticity')
axis equal tight
set(gca,'XTick',[],'YTick',[])

%%
% Plot the spectrum
disp('Plotting the spectrum');
figure('Name','Spectrum')
koopman.plotSpectrum( Spectrum, Amps )

%%
% Plotting the modes
figure('Name','Modes');
plotModes(4, Spectrum, Amps, Modes)

end

function plotModes( nPlotted, Spectrum, Amps, Modes )

% visualize only modes with positive angles -- to avoid replication
  sel = angle(Spectrum) >= 0;
  Spectrum = Spectrum(sel);
  Amps = Amps(sel);
  Modes = Modes(:,:,:,sel);

  SpeedModes = squeeze(Modes(1,:,:,:));
  VortModes = squeeze(Modes(2,:,:,:));

  nPlotted = 4;
  for k = 1:nPlotted
    subplot(nPlotted, 2, 2*k-1);
    pcolor( SpeedModes(:,:,k) );
    shading interp;
    title({sprintf('%d: Growth=%.1f, Freq=%.1f', k, real(Spectrum(k)),...
                   imag(Spectrum(k)) ), ...
           'Speed'});
    set(gca,'XTick',[],'YTick',[])
    axis equal tight

    subplot(nPlotted, 2, 2*k);
    pcolor( VortModes(:,:,k) );
    shading interp;
    title({sprintf('log_{10}(Amplitude)=%.1f, Phase=%.1f \\pi', log10(abs(Amps(k))),...
                   angle(Amps(k))/k ), ...
           'Vorticity'});
    set(gca,'XTick',[],'YTick',[])
    axis equal tight
  end

end