function [peakOmega, peakPower, pSpec] = getPowerSpectrum( y, dt )
%GETSPECTRUM Compute log-power spectrum of signal y and compute peaks.
%
% [PEAKOMEGA, PEAKPOWER, PSPEC] = GETPOWERSPECTRUM(Y,DT)
%  Compute the single-sided power spectrum of the signal Y, sampled at
%  uniform rate DT.
%  PEAKOMEGA - ang. frequencies of the most prominent peaks
%  PEAKPOWER - power at the most prominent peaks
%  PSPEC     - the entire power spectrum
%  OMEGA     - angular frequencies corresponding to the power spectrum
%
%  If no outputs are requested, the power spectrum is plotted and labeled.

% Copyright 2015 under BSD license (see LICENSE file).

% minimum prominence at which peak is detected (height to neighbor, in dB)
minProminence = 5;

% Sampling frequency
Fs = 1/dt;

% Compute DFT
NFFT = length(y);
Y = fft(y,NFFT);

% Compute power spectrum (norm of DFT)
pSpec = 10*log10( (abs(Y).^2)/NFFT );

% Truncate to single-sided DFT
pSpec = pSpec(1:round(NFFT/2));

% Compute physical frequency vector
F = ((0:1/NFFT:1-1/NFFT)*Fs);
% Convert to angular frequency
F = F(1:round(NFFT/2));
Omega = 2*pi*F;

% convert to row-vectors
pSpec = pSpec(:).';
Omega = Omega(:).';

% find peaks that are at least 10dB larger than neighbors
[peakPower,peakOmega] = findpeaks(pSpec,Omega,'SortStr','descend','NPeaks',5, ...
                        'MinPeakProminence',minProminence);

% Plot (using findpeaks internal plotting) and label peaks if there are no
% arguments requested
if nargout == 0

  findpeaks(pSpec,Omega,'SortStr','descend','NPeaks',5,'MinPeakProminence',minProminence)
  ylabel('dB')
  xlabel('Angular frequency')

  % place labels containing frequencies
  labels = arrayfun(@(loc)sprintf('%.1f',loc), ...
                    peakOmega,'UniformOutput',false);
  text(peakOmega+.02,peakPower,labels,'FontSize',7);

end
