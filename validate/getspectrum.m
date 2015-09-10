function [peak, Pyy] = getspectrum( y, dt )

Fs = 1/dt;
NFFT = length(y);
Y = fft(y,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs);
Pyy = Y.*conj(Y)/NFFT;
Om = 2*pi*F;

Pyy = Pyy(1:round(NFFT/2));
Om = Om(1:round(NFFT/2));
Pyy = Pyy(:).';
Om = Om(:).';

[~,sel] = findpeaks(Pyy);

sel = sel(:);

if nargout < 1
  semilogy(Om,Pyy,'x-')
  hold on
  semilogy( Om(sel), Pyy(sel), 'ro' );
  hold off;
  title('Power spectral density')
  xlabel('Angular frequency')
end

peak = sortrows( [Om(sel),Pyy(sel)], [-2,1] );
