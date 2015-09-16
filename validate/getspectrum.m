function [lsor, psor, Pyy] = getspectrum( y, dt )

% Sampling frequency
Fs = 1/dt;

% Compute DFT
NFFT = length(y);
Y = fft(y,NFFT);

% Compute power spectrum (norm of DFT)
Pyy = 10*log10( (abs(Y).^2)/NFFT );

% Truncate to single-sided DFT
Pyy = Pyy(1:round(NFFT/2));

% Compute physical frequency vector
F = ((0:1/NFFT:1-1/NFFT)*Fs);
% Convert to angular frequency
F = F(1:round(NFFT/2));
Om = 2*pi*F;

% convert to row-vectors
Pyy = Pyy(:).';
Om = Om(:).';

%
prominence = 5;

% find peaks that are at least 10dB larger than neighbors
[psor,lsor] = findpeaks(Pyy,Om,'SortStr','descend','NPeaks',5, ...
                        'MinPeakProminence',prominence);

% Plot and label
findpeaks(Pyy,Om,'SortStr','descend','NPeaks',5,'MinPeakProminence',prominence)

labels = arrayfun(@(loc)sprintf('%.1f',loc), ...
                  lsor,'UniformOutput',false);

text(lsor+.02,psor,labels,'FontSize',7);

ylabel('dB')
xlabel('Angular frequency')
